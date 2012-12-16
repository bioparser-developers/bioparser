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
 * The module parser.fasta is a set of function to parse fasta file format
 * See_Also:
 * Format specification http://en.wikipedia.org/wiki/FASTA_format
 * Copyright: Copyright Jonathan MERCIER  2012.
 *
 * License:   GPLv3+
 *
 * Authors:   Jonathan MERCIER aka bioinfornatics
 *
 * Source: parser/fasta.d
 */
 module parser.fasta;

import std.mmfile       : MmFile;
import std.array        : split, array, appender;
import std.algorithm    : filter, map, count, reduce;
import std.ascii        : isWhite, newline, toUpper;
import std.conv         : to;
import std.typecons     : Tuple;
import std.string       : stripRight;
import std.parallelism  : parallel;
import std.functional   : unaryFun;

import misc.ascii       : isNewline;

struct Fasta{
    string header;
    string[] sequence;

    this( string header, string[] sequence ){
        this.header     = header;
        this.sequence   = sequence;
    }

    @property string toString() const{
        string result = "";
        auto fasta  = appender(&result);
        fasta.put( header.idup ~ newline );
        foreach( seq; sequence )
            fasta.put( seq.idup ~ newline );
        return result;
    }

    bool opEquals( ref const Fasta o){
        return o.header == this.header && o.sequence == this.sequence;
    }
}

struct Section {
    private:
        MmFile _mappedFile;
        size_t _startHeader;
        size_t _endHeader; // endHeader + 1 should be start sequence then no eed to store it
        size_t _endSequence;

    public:
        this( MmFile mappedFile, size_t startHeader, size_t endHeader, size_t endSequence ){
            _mappedFile     = mappedFile;
            _startHeader    = startHeader;
            _endHeader      = endHeader;
            _endSequence    = endSequence;
        }

        @property string header() {
            debug writeln( "-> Section.header()" );
            return cast(string) _mappedFile[ _startHeader .. _endHeader ];
        }

        @property string[] sequence(){
            debug writeln( "-> Section.sequence()" );
            return (cast(string)_mappedFile[ _endHeader+1 .. _endSequence ])
                                                            .stripRight
                                                            .split(newline)
                                                            .map!( (s) => s.filter!( c => !isWhite( c ) )
                                                                            .map!( c => toUpper(c) ) // maybe is useless cost too cpu time
                                                                            .to!(string) )
                                                            .array
                                                            .to!(string[]);

        }

        Fasta toFasta() {
            return Fasta( header, sequence );
        }
}


/**
 * bySection is a way to iterate over a fasta file.
 *
 * Returns: A Section object with all information at each step
 *
 * Examples:
 * --------------------
 * import parser.fasta;
 * string filePath  = "~/myData.fastq";
 * foreach( count, section ; filePath.bySection() )
 *     Fasta fasta = section.toFasta();
 * --------------------
 */
struct bySection {
    private:
        MmFile      _mappedFile;
        size_t      _position;
        size_t      _number;
        Section     _current;
        bool        _isEmpty;

        void moveToNewline( ){
            debug writeln( "-> bySection.moveToNewline()" );
            bool isSearching = true;
            while( isSearching ){
                if( _position > _mappedFile.length || isNewline( cast(dchar)_mappedFile[_position] ) )
                    isSearching = false;
                else
                    _position++;
            }
            assert( isNewline( cast(dchar)_mappedFile[_position] ), "Fasta header found but never reach the header end." );
        }

        void moveToNextSection( ){
            debug writeln( "-> bySection.moveToNextSection()" );
            bool isSearching = true;
            while( isSearching ){
                if( _position >= _mappedFile.length || _mappedFile[_position] == 62 )
                    isSearching = false;
                else
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
                moveToNextSection();
                size_t      startHeader;
                size_t      endHeader; // endHeader + 1 should be start sequence then no need to store it
                size_t      endSequence;
                startHeader = _position;
                moveToNewline();
                endHeader   = _position;
                moveToNextSection();
                endSequence = _position;
                _current    = Section( _mappedFile, startHeader, endHeader, endSequence );
            }
        }
}


/**
 * doIndex is a way to index a fasta file. This allow to :
 *   - access to a random section as an array
 *   - iterate over the file in parallel
 *   - filter to element to use
 *
 * Returns: A Section object list
 *
 * Examples:
 * --------------------
 * import parser.fastqmyData.fasta";
 * Section[] sections  = doIndex!(s => s.identifier != "")( fastaFile ) // keep section where identifier is not empty
 *
 * Fasta fasta5 = sections[5].toFasta();
 *
 * foreach( index; parallel( iota(sections.length)  ) )
 *     Fasta fasta = sections[index].toFasta();
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
    string testFilename( string file = __FILE__, size_t line = __LINE__ ){
        return text("deleteme.", baseName(file), ".", line);
    }
}

unittest{
    writeln( "Testing parser.fasta module" );
    string fastaFile = testFilename() ;
    scope(exit) remove( fastaFile );
    scope(failure) writeln("module parser.fasta [NO]" );
    scope(success) writeln( "module parser.fasta [OK]" );
    write( fastaFile,
">SEQUENCE_1
MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG
LVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHK
IPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTL
MGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL
>SEQUENCE_2
SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI
ATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH
" );

    Fasta fasta1;
    Fasta fasta2;

    fasta1.header   = ">SEQUENCE_1";
    fasta1.sequence =   [
                            "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG",
                            "LVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHK",
                            "IPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTL",
                            "MGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL"
                        ];
    fasta2.header   = ">SEQUENCE_2";
    fasta2.sequence = ["SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI","ATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH"];

    Fasta[] fastaList = [ fasta1, fasta2 ];

    writeln( "\t+ Testing fasta file indexing" );
    Section[] sections = doIndex!( c => c.header != "" )( fastaFile );

    writeln( "\t+ Testing foreach loop over a fasta file" );
    foreach( i, section; fastaFile.bySection() ){
        Fasta fasta = section.toFasta();
        assert(fasta.header   == fastaList[i].header    , text( "Error fasta n°", i, " do not have the right header" ) );
        assert(fasta.sequence == fastaList[i].sequence  , text( "Error fasta n°", i, " do not have the right header" ) );
    }
    writeln("\t+ Testing to get fasta from his index" );
    assert( sections[0].toFasta().header   == fastaList[0].header    , "Error fasta n°0 do not have the right header" );
    assert( sections[0].toFasta().sequence == fastaList[0].sequence  , "Error fasta n°0 do not have the right sequence" );
    writeln("\t+ Testing to get fasta using retro range" );
    assert( sections.retro[0].toFasta().header   == fastaList.retro[0].header    , "Error fasta n°1 do not have the right header" );
    assert( sections.retro[0].toFasta().sequence == fastaList.retro[0].sequence  , "Error fasta n°1 do not have the right sequence" );
    writeln( "\t+ Testing parallel foreach loop over a fasta file" );
    foreach( i; parallel( iota( sections.length )  ) ){
        Fasta fasta = sections[i].toFasta();
        assert(fasta.header   == fastaList[i].header    , text( "Error fasta n°", i, " do not have the right header" ) );
        assert(fasta.sequence == fastaList[i].sequence  , text( "Error fasta n°", i, " do not have the right header" ) );
    }
}
