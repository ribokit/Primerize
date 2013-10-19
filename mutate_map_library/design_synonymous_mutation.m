function [syn_names, syn_seqs] = design_synonymous_mutation (sequence, seq_offset, cod_offset, std_dir, mode_flag, max_cod_mut, max_total_mut)

if nargin == 0; help( mfilename ); return; end;

if ~exist('seq_offset','var') || isempty(seq_offset); seq_offset = 0; end;
if ~exist('cod_offset','var') || isempty(cod_offset); cod_offset = 0; end;
if ~exist('std_dir','var') || isempty(std_dir) || (std_dir ~= 1 && std_dir ~= -1); std_dir = 1; end;

if ~exist('mode_flag','var') || isempty(mode_flag); mode_flag = {'123','only'}; end;
% read flag
mode_str = lower([mode_flag{1:end}]);
multi_flag = [];
if strfind(mode_str, '1'); multi_flag = [multi_flag, 1]; end;
if strfind(mode_str, '2'); multi_flag = [multi_flag, 2]; end;
if strfind(mode_str, '3'); multi_flag = [multi_flag, 3]; end;
if strfind(mode_str, 'full');
    full_flag = 1;
else
    full_flag = 0;
end;

if ~exist('max_cod_mut','var') || isempty(max_cod_mut); max_cod_mut = 4; end;
if ~isinteger(max_cod_mut); max_cod_mut = round(max_cod_mut); end;
if sign(max_cod_mut) == -1; max_cod_mut = 0; end;
if ~exist('max_total_mut','var') || isempty(max_total_mut); max_total_mut = 0; end;
if ~isinteger(max_total_mut); max_total_mut = round(max_total_mut); end;
if sign(max_total_mut) == -1; max_total_mut = 0; end;

% parse all codons
cod_tab = codon_table;
[~, codons] = show_ORF(sequence, cod_offset, std_dir);

syn_names = {};
syn_seqs = {};
syn_codons = {};
% index of output mutations cell
mut_count = 0;
for i = 1:length(codons)
    
    % get list of target mutation codons
    codon_mutations = find_codon_from_codon(codons{i}, cod_tab);
    for j = 1:length(codon_mutations)
        name_count = 0;
        mut = {};
        % make the mutation set cell of names
        for k = 1:3
            if ~strcmp(codon_mutations{j}{1}(k),codons{i}(k));
                name_count = name_count + 1;
                % numbering for strand direction
                if std_dir == 1
                    mut{name_count} = [codons{i}(k),...
                        num2str(seq_offset + cod_offset + (i-1)*3 + k),...
                        codon_mutations{j}{1}(k)];
                else
                    mut{name_count} = [strrep(complement(codons{i}(k)),'T','U'),...
                        num2str(seq_offset + length(sequence) - cod_offset - ((i-1)*3 + k -1)),...
                        strrep(complement(codon_mutations{j}{1}(k)),'T','U')];
                end;
            end;
        end;
        
        % screen only the mutants specified by mode_flag
        if ismember(length(mut), multi_flag);
            mut_count = mut_count + 1;
            syn_names(mut_count) = {mut};
            syn_seqs(mut_count) = design_mutants({}, sequence, syn_names(mut_count), -seq_offset, '','',0);
            for k = 1:length(codons)
                if k == i;
                    syn_codons(mut_count,k) = codon_mutations{j};
                else
                    syn_codons(mut_count,k) = codons(k);
                end;
            end;
        end;
    end;
end;

% output all mutants
fprintf('Output\n======\n');
% add WT to first row because output_triplet_seq requires that
syn_seqs_temp = [sequence, syn_seqs];
syn_codons_temp = cell(size(syn_codons,1)+1,length(sequence(1+cod_offset: length(sequence)- mod((length(sequence)- cod_offset),3)))/3);
syn_codons = [rot90(codons); syn_codons];
for i = 1:size(syn_codons,1)
    for j = 1:size(syn_codons,2)
        syn_codons_temp(i,j) = syn_codons(i,j);
    end;
end;
output_triplet_seq(syn_names, syn_codons_temp, syn_seqs_temp, cod_offset, 0, std_dir);
fprintf(['Total = ', num2str(length(syn_names)),'; excluding WT.\n']);

