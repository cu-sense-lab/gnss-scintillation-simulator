function DateStamp=generateDateStamp
%Unique output file name for setup parameters and PropCodeSim output
    vec = fix(clock);
    hours   = num2str( vec(4) );
    minutes = num2str( vec(5) );
    seconds = num2str( vec(6) );
    DateStamp = ['_' date '_H' hours '_M' minutes '_S' seconds];