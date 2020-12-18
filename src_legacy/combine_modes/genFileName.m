function name = genFileName( fh, index )

    markPos = strfind( fh.name, '(*)' );
    if fh.padwidth == 0
	    name = strcat( fh.name(1:markPos-1), ...
	        num2str( index ), fh.name(markPos+3:end) );
	else
	    name = strcat( fh.name(1:markPos-1), ...
	        num2str( index, strcat( '%0.', num2str(fh.padwidth), 'd' ) ), ...
	        fh.name(markPos+3:end) );
	end

end