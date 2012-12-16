module misc.ascii;

bool isNewline(dchar c) @safe pure nothrow {
    return ( c == 0x0A || c == 0x0D );
}
version(unittest){
    import std.stdio         : writeln;
}

unittest{
    writeln( "Testing misc.ascii module" );
    writeln( "\t+ Testing \\n character" );
    assert( '\n'.isNewline  , "Error function do not recognize \\n character"       );
    writeln( "\t+ Testing \\r character" );
    assert( '\r'.isNewline  , "Error function do not recognize \\r character"       );
    writeln( "\t+ Testing space character" );
    assert( ! ' '.isNewline , "Error function should not recongize space character" );
}
