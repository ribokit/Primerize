function primer_sequences = output_primers_server( primers_all, sequence, tag )

if ~exist( 'tag','var' ); tag = 'primer'; end;
OUTPUT_STAGGER = 1;

% Assume primers are alternating in directionality.
% count = 0;

num_primers = size( primers_all,2);

blank_line = '';
COLWIDTH = 138;
for k = 1:max(length(sequence),COLWIDTH)
    blank_line(k) = ' ';
end;

seq_line_prev = blank_line;
bp_lines = cell(1, num_primers);
seq_lines = cell(1, num_primers);
primer_sequences = cell(1, num_primers);

for j = 1:num_primers
    
    primers = primers_all( :, j );
    seq_line = blank_line;
    
    seg_start = primers(1);
    seg_end = primers(2);
    direction = primers(3);
    
    if ( direction == 1 );
        for k = seg_start:seg_end;
            seq_line(k) = sequence(k);
        end;
        primer_sequences{j} = sequence( seg_start:seg_end );
        
        if ( seg_end + 1 <= length(sequence ));
            seq_line( seg_end+1 ) = '-';
        end;
        if ( seg_end + 2 <= length(sequence ));
            seq_line( seg_end+2 ) = '>';
        end;
        
        text_out = num2str( j ) ;
        if ( seg_end + 2 + length(text_out)  <= length( sequence) );
            seq_line( seg_end + 2 + [1:length(text_out)] ) = text_out;
        end;
    else
        for k = seg_start:seg_end;
            seq_line(k) = reverse_complement( sequence(k) );
        end
        primer_sequences{j} = reverse_complement(sequence( seg_start:seg_end ));
        
        if ( seg_start - 1 >= 1 );
            seq_line( seg_start-1 ) = '-';
        end;
        if ( seg_start - 2 >= 1 );
            seq_line( seg_start-2 ) = '<';
        end;
        
        text_out = num2str( j ) ;
        if ( seg_start - 2 - length(text_out)  >= 1 )
            seq_line( seg_start - 3 - length(text_out) + [1:length(text_out)] ) = text_out;
        end;
    end;
    
    bp_line = blank_line;
    overlap_seq = '';
    last_bp_pos = 1;
    for k = 1:length(sequence)
        if ( isempty(strfind( 'ACGT',seq_line_prev(k))) || ...
                isempty(strfind( 'ACGT',seq_line(k))) );
            bp_line(k) = ' ';
        else
            bp_line(k) = '|';
            last_bp_pos = k;
            overlap_seq = [ overlap_seq, sequence(k) ];
        end;
    end;
    
    if (last_bp_pos > 1);
        Tm = calc_Tm( overlap_seq, 0.2e-6,0.1,0.0015);
        Tm_text = num2str( Tm, '{%2.1f}' );
        bp_line( last_bp_pos + [1:length(Tm_text)+1] ) = [' ',Tm_text];
    end;
    bp_lines{j} = bp_line;
    seq_lines{j} = seq_line;
    %fprintf( '%s\n%s\n', bp_line, seq_line );
    seq_line_prev = seq_line;
end;

fprintf('#\n');
fprintf('primers\t\tlength\tsequences\n');
for j = 1:length( primer_sequences )
    fprintf( '%s-%d\t%d\t%s\n', tag, j, length(primer_sequences{j}), primer_sequences{j} );
end;
fprintf('#\n');

if ~OUTPUT_STAGGER;
    fprintf('%s',sequence);
    for k = 1:length(seq_lines)
        fprintf( '%s\n%s\n',bp_lines{k}, seq_lines{k} );
    end;
    fprintf('%s\n\n',complement(sequence));
else
%     blank_line = blank_line( 1:COLWIDTH);
    for n = 1: floor( (length(sequence)-1)/ COLWIDTH)+1
        start_pos = COLWIDTH*(n-1) + 1;
        end_pos   = min( COLWIDTH*(n-1) + COLWIDTH, length(sequence));
        out_line = sequence(start_pos:end_pos);
        fprintf( '~%s\n', out_line);
        for k = 1:length(seq_lines)
            if ~isempty(strrep(bp_lines{k}(end_pos:end),' ','')) && isempty(strfind(strrep(bp_lines{k}(end_pos:end),' ',''),'|')) && isempty(strrep(bp_lines{k}(1:start_pos-1),' ',''));
                bp_line = deblank(bp_lines{k} (start_pos:end ));
            elseif isempty(strfind(bp_lines{k}(start_pos:end_pos),'|'));
                bp_line = repmat(' ',1,end_pos-start_pos+1);
            else
                bp_line = bp_lines{k} (start_pos:end_pos );
            end;                
            seq_line = seq_lines{k}(start_pos:end_pos );
            if ~isempty(strrep(bp_line,' ','')) || ~isempty(strrep(seq_line,' ',''));
                fprintf( '$%s\n', bp_line);
                if mod(k,2) == 0;
                    fprintf('!');
                else
                    fprintf('^');
                end;
                fprintf('%s\n', seq_line);
            end;
        end;
        fprintf( '$%s\n=%s\n\n\n', repmat(' ', 1, end_pos - start_pos + 1), complement(out_line));
    end;
end;


