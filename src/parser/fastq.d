/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 */

/**
 * The module parser.fastq is a set of function to parse fastq file format
 * See_Also:
 * http://en.wikipedia.org/wiki/FASTQ_format
 * Copyright: Copyright Jonathan MERCIER  2012.
 *
 * License:   GPLv3+
 *
 * Authors:   Jonathan MERCIER aka bioinfornatics
 *
 * Source: parser/fastq.d
 */
 module parser.fastq;

import std.mmfile       : MmFile;
import std.array        : split, array, appender;
import std.algorithm    : filter, map, count, reduce;
import std.ascii        : isWhite, newline, toUpper;
import std.conv         : to;
import std.string       : stripRight, format;
import std.typecons     : Tuple;
import misc.ascii       : isNewline;
import std.functional   : unaryFun;
debug import std.stdio;

struct Fastq {
    string identifier;
    string sequence;
    string description;
    string quality;

    this( string identifier, string sequence, string description, string quality ){
        this.identifier     = identifier;
        this.sequence       = sequence;
        this.description    = description;
        this.quality        = quality;
    }

    @property string toString() const{
        string result;
        auto fastq  = appender(&result);

        fastq.put( "%s%s".format( identifier, newline ) );

        fastq.put(  sequence );
        fastq.put( newline );
        fastq.put( "%s%s".format( description, newline ) );

        fastq.put( quality );

        return result;
    }

    bool opEquals( ref const Fastq o){
        return  o.identifier    == identifier   &&
                o.sequence      == sequence     &&
                o.description   == description  &&
                o.quality       == quality;
    }
}


struct Section {
    private:
        MmFile _mappedFile;
        size_t _startIdentifier;
        size_t _startSequence;
        size_t _startDescription;
        size_t _startQuality;
        size_t _endQuality;

    public:
        this(MmFile mappedFile, size_t startIdentifier, size_t startSequence, size_t startDescription, size_t startQuality, size_t endQuality ){
            _mappedFile      = mappedFile;
            _startIdentifier = startIdentifier;
            _startSequence   = startSequence;
            _startDescription= startDescription;
            _startQuality    = startQuality;
            _endQuality      = endQuality;
        }

        @property string identfier(){
            debug writeln( "-> Section.identfier()" );
            return cast(string) _mappedFile[ _startIdentifier .. _startSequence - 1 ]; // if a problem occure with _r\ switch by using stripRingth instead -1
        }

        @property string sequence(){
            debug writeln( "-> Section.sequence()" );
            return (cast(string) _mappedFile[ _startSequence .. _startDescription - 1 ])
                                                    .filter!( c => !isWhite( c ) )
                                                    .map!( c => toUpper(c) )
                                                    .to!(string);
        }

        @property string description(){
            debug writeln( "-> Section.description()" );
            return cast(string) _mappedFile[ _startDescription .. _startQuality - 1 ];
        }

        @property string quality(){
            debug writeln( "-> Section.quality()" );
            return (cast(string) _mappedFile[ _startQuality .. _endQuality ])
                                                    .filter!( c => !isWhite( c ) )
                                                    .map!( c => toUpper(c) )
                                                    .to!(string);
        }

        Fastq toFastq(){
            debug writeln( "-> Section.toFastq()" );
            return Fastq(identfier, sequence, description, quality );
        }
}


/**
 * bySection is a way to iterate over a fastq file.
 *
 * Returns: A Section object with all information at each step
 *
 * Examples:
 * --------------------
 * import parser.fastq;
 * string filePath  = "~/myData.fastq";
 * foreach( count, section ; filePath.bySection() )
 *     Fastq fastq = section.toFastq();
 * --------------------
 */
struct bySection {
    private:
        MmFile      _mappedFile;
        size_t      _position;
        size_t      _number;
        size_t      _letters;
        Section     _current;
        bool        _isEmpty;

        void moveToNextSection(){
            debug writeln( "-> bySection.moveToNextSection()" );
            while( _position < _mappedFile.length &&  _mappedFile[_position] != 64 )
                _position++;
        }

        void moveToEndSequence(){
            debug writeln( "-> bySection.moveToEndSequence()" );
            while( _mappedFile[_position] != 43) { // 43 ==> +
                assert( _position < _mappedFile.length, "Error index out of file" );
                if( ! isWhite( cast(dchar) _mappedFile[_position] ) )
                    _letters++;
                _position++;
            }
        }

        void moveToNewline(){
            debug writeln( "-> bySection.moveToNewline()" );
            while( ! isNewline( cast(dchar)_mappedFile[_position] ) )
                _position++;
        }

        void moveToEndQuality(){
            debug writeln( "-> bySection.moveToEndQuality()" );
            while( _letters  != 0 ){
                assert( _position < _mappedFile.length, "Error index out of file" );
                if( ! isWhite( cast(dchar) _mappedFile[_position] ) )
                    _letters--;
                _position++;
            }
        }

    public:
        this( in string filename ){
            _mappedFile = new MmFile( filename );
            if( _mappedFile.length != 0 )
                popFront();
        }

        @property bool empty() const{
            debug writeln( "-> bySection.empty()" );
            return _isEmpty;
        }

        @property size_t length() const {
            debug writeln( "-> bySection.length()" );
            return _mappedFile.length;
        }

        alias length opDollar;

        @property Tuple!( size_t, Section ) front() {
            debug writeln( "-> bySection.front()" );
            Tuple!( size_t, Section ) result;
            result[0] = _number;
            result[1] = _current;
            _number++;
            return result;
        }

        void popFront() {
            debug writeln( "-> bySection.popFront()" );
            _isEmpty = ( _mappedFile.length == 0 || _mappedFile.length == _position );
            if( ! _isEmpty ){
                size_t startIdentifier;
                size_t startSequence;
                size_t startDescription;
                size_t startQuality;
                size_t endQuality;
                moveToNextSection();
                startIdentifier = _position;
                moveToNewline(); // end identifier
                startSequence = _position + 1;
                debug writeln( cast(string) _mappedFile[ startIdentifier .. startSequence ] );
                moveToEndSequence();
                startDescription = _position + 1;
                debug writeln( cast(string) _mappedFile[ startSequence .. startDescription ] );
                moveToNewline(); // end description
                startQuality = _position + 1;
                debug writeln( cast(string) _mappedFile[ startDescription .. startQuality ] );
                moveToEndQuality();
                endQuality = _position;
                debug writeln( cast(string) _mappedFile[ startQuality .. endQuality ] );
                _current = Section( _mappedFile, startIdentifier, startSequence, startDescription, startQuality, endQuality );
            }
        }
}


/**
 * doIndex is a way to index a fastq file. This allow to :
 *   - access to a random section as an array
 *   - iterate over the file in parallel
 *   - filter to element to use
 *
 * Returns: A Section object list
 *
 * Examples:
 * --------------------
 * import parser.fastqmyData.fastq";
 * Section[] sections  = doIndex!(s => s.identifier != "")( fastqFile ) // keep section where identifier is not empty
 *
 * Fastq fastq5 = sections[5].toFastq();
 *
 * foreach( index; parallel( iota(sections.length)  ) )
 *     Fastq fastq = sections[index].toFastq();
 * --------------------
 */
template doIndex(alias pred) if ( is( typeof( unaryFun!pred ) ) ) {
    Section[] doIndex( in string filename ){
        Section[]   sections;
        auto        list = appender( &sections );

        foreach( i, section; filename.bySection() ){
            if( unaryFun!pred( section ) )
                list.put( section );
        }

        return sections;
    }
}


version(unittest){
    import std.range        : retro, zip, iota;
    import std.stdio        : writeln;
    import std.file         : remove, write;
    import std.path         : baseName;
    import std.conv         : text;
    import std.parallelism  : parallel;
    string testFilename(string file = __FILE__, size_t line = __LINE__){
        return text("deleteme.", baseName(file), ".", line);
    }
}

unittest {
    writeln( "Testing parser.fastq module" );
    string fastqFile = testFilename();
    scope(exit) remove( fastqFile );
    scope(failure) writeln("module parser.fasta [NO]" );
    scope(success) writeln( "module parser.fasta [OK]" );
    write( fastqFile,
"@SOUFRE:176:C0BULACXX:6:1101:1109:2127/1
TG
TTTTCAAAGTCGTCGCCCTGACGAGTTGTCCCCAACATTGTGGTCCACGGAATGAGGTTCATGCGGGGGCTATTTCCGTCTCGACTTGTCAACTTGAAA
+ testing
B@BFFF FFHFHFHJIJJJJJJJGIIIJJIIJJJG
JIJJJJJJIIIFHIJJJJIGHEHHHCDFFEEEDDDDDDBBDDD:@?CDDDD@@BBDDDEEDDDC@>:
@SOUFRE:176:C0BULACXX:6:1101:1701:2131/1
CGTCGAAATATTCGCGGGGCTCGCCAGTGTTGTCGAACACCACTCTGCCGTCTTGTTCGCGAATCAAGCGCTTGTCGAGTACGAGTTTCAAGCATGCGGCG
+
CCCFFFFFGHHHHJJJJJJJJJIJJGHEFFEDFEDDDBBCD<BDDDDDDD@DDDCDDCDDDDBDDDDDDDDDDDDDD?B>?CDBB2?CDEDDDDCDDDDBD
@SOUFRE:176:C0BULACXX:6:1101:1579:2169/1
CACCGCTCGGTGCGCCGCTGTGCGCCGCCATCCCGCGCGACGGTGCCGGC
+
@??BD?DFDHC;?EAFGDGIFHIHID;FGBCE>E=6=AB=;@B:<@CBBB"
    );

    Fastq fastq1;
    Fastq fastq2;
    Fastq fastq3;

   fastq1.identifier    = "@SOUFRE:176:C0BULACXX:6:1101:1109:2127/1";
   fastq1.sequence      = "TGTTTTCAAAGTCGTCGCCCTGACGAGTTGTCCCCAACATTGTGGTCCACGGAATGAGGTTCATGCGGGGGCTATTTCCGTCTCGACTTGTCAACTTGAAA";
   fastq1.description   = "+ testing";
   fastq1.quality       = "B@BFFF FFHFHFHJIJJJJJJJGIIIJJIIJJJGJIJJJJJJIIIFHIJJJJIGHEHHHCDFFEEEDDDDDDBBDDD:@?CDDDD@@BBDDDEEDDDC@>:";

   fastq2.identifier    = "@SOUFRE:176:C0BULACXX:6:1101:1701:2131/1";
   fastq2.sequence      = "CGTCGAAATATTCGCGGGGCTCGCCAGTGTTGTCGAACACCACTCTGCCGTCTTGTTCGCGAATCAAGCGCTTGTCGAGTACGAGTTTCAAGCATGCGGCG";
   fastq2.description   = "+";
   fastq2.quality       = "CCCFFFFFGHHHHJJJJJJJJJIJJGHEFFEDFEDDDBBCD<BDDDDDDD@DDDCDDCDDDDBDDDDDDDDDDDDDD?B>?CDBB2?CDEDDDDCDDDDBD";

   fastq3.identifier    = "@SOUFRE:176:C0BULACXX:6:1101:1579:2169/1";
   fastq3.sequence      = "CACCGCTCGGTGCGCCGCTGTGCGCCGCCATCCCGCGCGACGGTGCCGGC";
   fastq3.description   = "+";
   fastq3.quality       = "@??BD?DFDHC;?EAFGDGIFHIHID;FGBCE>E=6=AB=;@B:<@CBBB";

    Fastq[] fastqList   = [ fastq1, fastq2, fastq3 ];
    writeln( "\t+ Testing fastq file indexing" );
    Section[] sections  = doIndex!( c => c.identfier != "" )( fastqFile );

    writeln( "\t+ Testing foreach loop over a fastq file" );
    foreach( index, section; fastqFile.bySection() ){
        Fastq fastq = section.toFastq();
        assert( fastq.identifier == fastqList[index].identifier  , text( "Error fastq n°", index, " do not have the right identifier") );
        assert( fastq.sequence   == fastqList[index].sequence    , text( "Error fastq n°", index, " do not have the right sequence"  ) );
    }
    writeln("\t+ Testing to get fastq from his index" );
    assert( sections[0].toFastq().identifier   == fastqList[0].identifier  , "Error fastq n°0 do not have the right identifier");
    assert( sections[0].toFastq().sequence     == fastqList[0].sequence    , "Error fastq n°0 do not have the right sequence"  );
    writeln("\t+ Testing to get fastq using retro range" );
    assert( sections.retro[0].toFastq().identifier == fastqList.retro[0].identifier  , "Error fastq n°0 do not have the right identifier");
    assert( sections.retro[0].toFastq().sequence   == fastqList.retro[0].sequence    , "Error fastq n°0 do not have the right sequence"  );
    writeln( "\t+ Testing parallel foreach loop over a fastq file" );
    foreach( index; parallel( iota(sections.length)  ) ){
        Fastq fastq = sections[index].toFastq();
        assert(fastq.identifier  == fastqList[index].identifier  , text( "Error fastq n°", index, " do not have the right identifier") );
        assert(fastq.sequence    == fastqList[index].sequence    , text( "Error fastq n°", index, " do not have the right sequence"  ) );
    }
}
