module misc.stats;

import std.typecons     : Tuple;
import std.parallelism  : parallel;

alias Tuple!(size_t[char], "quantity", size_t, "total") Result;

@safe pure
Result lettersCounter( const ref string sequence ){
    Result result;
    size_t length = 0;
    foreach( letter; sequence ){
        if( letter in result[0]) result[0][letter]++;
        else result[0][letter] = 1;
    }
    result[1] = sequence.length;
    return result;
}

Result lettersCounter( const ref string[] sequences ){
    Result   result;
    Result[] baseNumbers  = new Result[](sequences.length);
    foreach(index, ref seq; parallel(sequences))
        baseNumbers[index] = lettersCounter( seq );
    foreach(baseNumber; baseNumbers){
        foreach( key, value; baseNumber[0] ){
            if( key in result[0] ) result[0][key] += value;
            else result[0][key] = value;
        }
        result[1] += baseNumber[1];
    }
    return result;
}

version(unittest){
    import std.range;
    import std.stdio         : writeln;
    import std.file          : remove, write;
    import std.path          : baseName;
    import std.conv          : text;
    string testFilename(string file = __FILE__, size_t line = __LINE__){
        return text("deleteme.", baseName(file), ".", line);
    }
}
unittest{
    writeln( "Testing misc.stats module" );
    string      sequence1   = "ACGATCGATCGATCGATCGATCGATCGATCGATCGGAGATAG";
    string[]    sequenceList= ["AGCTAGCTGATCGATCGTAGCTGATGCA","CCGTAGCTAGCTATTAGCTAGCCTCGTAGCTAGCTAGCTA","CGATCGATCGATCGATGCTATGCG"];

    writeln( "\t+ Testing  lettersCounter over a string sequence" );
    Result r1 = sequence1.lettersCounter;
    assert( r1[0]['A']      == 12   , text("Error found ", r1[0]['A'], " instead of 12") );
    assert( r1.quantity['C']== 9    , text("Error found ", r1[0]['C'], " instead of 12") );
    assert( r1[0]['T']      == 9    , text("Error found ", r1[0]['G'], " instead of 9" ) );
    assert( r1[0]['G']      == 12   , text("Error found ", r1[0]['T'], " instead of 12") );
    assert( r1.total        == 42   , text("Error found ", r1[1]     , " instead of 42") );
    
    writeln( "\t+ Testing  lettersCounter over a list of string sequences" );
    Result r2 = sequenceList.lettersCounter;
    assert( r2[0]['A']      == 21   , text("Error found ", r2[0]['A'], " instead of 21") );
    assert( r2.quantity['C']== 23   , text("Error found ", r2[0]['C'], " instead of 23") );
    assert( r2[0]['T']      == 24   , text("Error found ", r2[0]['G'], " instead of 24") );
    assert( r2[0]['G']      == 24   , text("Error found ", r2[0]['T'], " instead of 24") );
    assert( r2.total        == 92   , text("Error found ", r2[1]     , " instead of 92") );

}
