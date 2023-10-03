clear;
tic
[user,sys] = memory;%place in the begging of your program

%        Species Description      GenBank Accession
data = {'Bacillus maritimus', 'KP317497';
    'Bacillus wakoensis', 'NR_040849';
    'Bacillus australimaris', 'NR_148787';
    'Bacillus xiamenensis', 'NR_148244';
    'Escherichia coli', 'J01859';
    'Streptococcus himalayensis', 'NR_156072';
    'Streptococcus halotolerans', 'NR_152063';
    'Streptococcus tangierensis', 'NR_134818';
    'Streptococcus cameli', 'NR_134817';
    'Thermus amyloliquefaciens', 'NR_136784';
    'Thermus tengchongensis', 'NR_132306';
    'Thermus thermophilus', 'NR_037066';
    'Thermus filiformis', 'NR_117152'
    %'Human' 'V00662';
    %         'Pygmy chimpanzee' 'D38116';
    %         'Common chimpanzee' 'D38113';
    %         'Gorilla' 'D38114';
    %         'Orangutan' 'D38115';
    %         'Gibbon' 'X99256';
    %         'Baboon' 'Y18001';
    %         'Horse' 'X79547';
    %         'White rhinoceros' 'Y07726';
    %         'Harbor seal' 'X63726';
    %         'Gray seal' 'X72004';
    %         'Cat' 'U20753';
    %         'Fin whale' 'X61145';
    %         'Blue whale' 'X72204';
    %         'Cow' 'V00654';
    %         'Rat' 'X14848';
    %         'Mouse' 'V00711';
    %         'Platypus' 'X83427';
    %'German_Neanderthal'      'D38116';
    %         'Russian_Neanderthal'     'AF254446';
    %         'European_Human'          'X90314'  ;
    %         'Mountain_Gorilla_Rwanda' 'AF089820';
    %         'Chimp_Troglodytes'       'AF176766';
    %         'Puti_Orangutan'          'AF451972';
    %         'Jari_Orangutan'          'AF451964';
    %         'Western_Lowland_Gorilla' 'AY079510';
    %         'Eastern_Lowland_Gorilla' 'AF050738';
    %         'Chimp_Schweinfurthii'    'AF176722';
    %         'Chimp_Vellerosus'        'AF315498';
    %         'Chimp_Verus'             'AF176731';
    };


%data E.Coli
data1 = {  'CB9615' 'CB9615';
    'APEC01' 'NZ_KB642078';
    'ATCC8739' 'ATCC8739';
    'B18BS512' 'B18BS512';
    'B4Sb227' 'B4Sb227';
    'BW2952' 'BW2952';
    'CFT073' 'CFT073';
    'D1Sd197' 'D1Sd197';
    'DH10B' 'DH10B';
    'E234869' 'E234869';
    'E24377A' 'E24377A';
    'ED1a' 'ED1a';
    'EDL933' 'EDL933';
    'F2a2457T' 'F2a2457T';
    'F2a301' 'F2a301';
    'F5b8401' 'F5b8401';
    'HS' 'HS';
    'IAI1' 'IAI1';
    'IAI39' 'IAI39';
    'MG1655' 'MG1655';
    'S88'  'S88';
    'Sakai' 'Sakai';
    'SE11' 'SE11';
    'SMS35' 'SMS35';
    'SSSs046' 'SSSs046';
    'UMN026' 'UMN026';
    'UTI89' 'UTI89'
    'W3110' 'W3110';
    
    };


%reading data from local benchmark databse subfolder
%dn='assembled-fish_mito\';
%dn='..\assembled-ecoli\';
%dn='..\unsimulated-ecoli_shigella\';
%dn='unsimulated-yersinia\';
%dn='16sRiboDNA\';
%dn='18EutherianMammal\';
%dn='21 HIV-1\';
dn='48 HEV\';
%dn='NADH\';

data=dir(strcat(dn,'*.fasta'));


descriptor=0;
descriptorCgrReduced=0;
prev=0;
partWiseDesc=0;


M=0;
surfTrack=1;
surfCount=1;


%dynamic k selection
lengthSequences=0;

for seqIter=1:length(data)
    [header, sequence]= fastaread(strcat(dn,data(seqIter).name));
    lengthSequences(seqIter)=length(sequence);
end

avgLength=mean(lengthSequences);
avgLength = int64(avgLength);

k=0;

if avgLength>=1 && avgLength<=20000
    k=4;
else
    k=6;
end

%k=10; 

for seqIter=1:length(data)
    % for getting data from gene bank
    %     primates(seqIter).Header   = data{seqIter,1};
    %     primates(seqIter).Sequence = getgenbank(data{seqIter,2},'sequenceonly','true');
    %     s= getgenbank(data{seqIter,2});
    %     sequence=primates(seqIter).Sequence;
    
    
    %for reading from local database subfolder
    %[header, sequence]= fastaread('D38116.fna');
    [header, sequence]= fastaread(strcat(dn,data(seqIter).name));
    
    %this two line only for 16s ribosomal DNA to ommit parts after space
%     newStr = split(header," ");
%     header=char(newStr{1});
%     header=[char(newStr{1}) ' ' char(newStr{2}) ' ' char(newStr{3})];
%     
    primates(seqIter).Header=header;
    primates(seqIter).Sequence=sequence;
    
    kmer=k; %take three length string
    %featureVector=CGRImageKmerCount(sequence,kmer);
    featureVector=CGRTwoDList(sequence,kmer);
    [row col]=size(featureVector);
    descriptor(seqIter,1:col)=featureVector;
    
    %disp(seqIter);
end


distances=pdist(descriptor);
UPGMAtree = seqlinkage(distances,'UPGMA',primates);

toc
[user2,sys2] = memory;%place in the end of your program
memory_used_in_bytes=user2.MemAvailableAllArrays-user.MemAvailableAllArrays;
disp(memory_used_in_bytes/(1024*1024));
h = plot(UPGMAtree,'orient','left');
title('UPGMA Distance Tree of Primates');
ylabel('Evolutionary distance');
%get canonical to identify the similarity
[Pointers, Distances, Names] = getcanonical(UPGMAtree);
%Create Newick-formatted character vector for benchmark
nwk = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk);
newNwkVect(1,1:nwkc)=nwk;


UPGMAtree = seqlinkage(distances,'single',primates);
h = plot(UPGMAtree,'orient','left');
title('single Distance Tree of Primates');
ylabel('Evolutionary distance');
%get canonical to identify the similarity
[Pointers, Distances, Names] = getcanonical(UPGMAtree);
%Create Newick-formatted character vector for benchmark
nwk1 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk1);
newNwkVect(2,1:nwkc)=nwk1;

UPGMAtree = seqlinkage(distances,'complete',primates);
h = plot(UPGMAtree,'orient','left');
title('complete Distance Tree of Primates');
ylabel('Evolutionary distance');
%get canonical to identify the similarity
[Pointers, Distances, Names] = getcanonical(UPGMAtree);
%Create Newick-formatted character vector for benchmark
nwk2 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk2);
newNwkVect(3,1:nwkc)=nwk2;

UPGMAtree = seqlinkage(distances,'weighted',primates);
h = plot(UPGMAtree,'orient','left');
title('weighted Distance Tree of Primates');
ylabel('Evolutionary distance');
%get canonical to identify the similarity
[Pointers, Distances, Names] = getcanonical(UPGMAtree);
%Create Newick-formatted character vector for benchmark
nwk3 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk3);
newNwkVect(4,1:nwkc)=nwk3;

UPGMAtree = seqlinkage(distances,'centroid',primates);
h = plot(UPGMAtree,'orient','left');
title('centroid Distance Tree of Primates');
ylabel('Evolutionary distance');
%get canonical to identify the similarity
[Pointers, Distances, Names] = getcanonical(UPGMAtree);
%Create Newick-formatted character vector for benchmark
nwk4 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk4);
newNwkVect(5,1:nwkc)=nwk4;

UPGMAtree = seqlinkage(distances,'median',primates);
% h = plot(UPGMAtree,'orient','left');
% title('median Distance Tree of Primates');
% ylabel('Evolutionary distance');
%get canonical to identify the similarity
[Pointers, Distances, Names] = getcanonical(UPGMAtree);
%Create Newick-formatted character vector for benchmark
nwk5 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk5);
newNwkVect(6,1:nwkc)=nwk5;

UPGMAtree1 = seqneighjoin(distances,'equivar',primates)
% h = plot(UPGMAtree1,'orient','left');
% title('UPGMA Distance Tree of Primates');
% ylabel('Evolutionary distance');
%get canonical to identify the similarity
[Pointers, Distances, Names] = getcanonical(UPGMAtree1);
%Create Newick-formatted character vector for benchmark
nwk6 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk6);
newNwkVect(7,1:nwkc)=nwk6;

UPGMAtree1 = seqneighjoin(distances,'firstorder',primates)
% h = plot(UPGMAtree1,'orient','left');
% title('UPGMA Distance Tree of Primates');
% ylabel('Evolutionary distance');
%get canonical to identify the similarity
[Pointers, Distances, Names] = getcanonical(UPGMAtree1);
%Create Newick-formatted character vector for benchmark
nwk7 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk7);
newNwkVect(8,1:nwkc)=nwk7;




%exp 2 dist squaredeuclidean + upgma
distances=pdist(descriptor,'squaredeuclidean');
UPGMAtree = seqlinkage(distances,'UPGMA',primates);
%Create Newick-formatted character vector for benchmark
nwk8 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk8);
newNwkVect(9,1:nwkc)=nwk8;

UPGMAtree = seqlinkage(distances,'single',primates);
%Create Newick-formatted character vector for benchmark
nwk9 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk9);
newNwkVect(10,1:nwkc)=nwk9;

UPGMAtree = seqlinkage(distances,'complete',primates);
%Create Newick-formatted character vector for benchmark
nwk10 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk10);
newNwkVect(11,1:nwkc)=nwk10;

UPGMAtree = seqlinkage(distances,'weighted',primates);
%Create Newick-formatted character vector for benchmark
nwk11 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk11);
newNwkVect(12,1:nwkc)=nwk11;

UPGMAtree = seqlinkage(distances,'centroid',primates);
%Create Newick-formatted character vector for benchmark
nwk12 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk12);
newNwkVect(13,1:nwkc)=nwk12;


UPGMAtree = seqlinkage(distances,'median',primates);
%Create Newick-formatted character vector for benchmark
nwk13 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk13);
newNwkVect(14,1:nwkc)=nwk13;

UPGMAtree1 = seqneighjoin(distances,'equivar',primates)
%Create Newick-formatted character vector for benchmark
nwk14 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk14);
newNwkVect(15,1:nwkc)=nwk14;

UPGMAtree1 = seqneighjoin(distances,'firstorder',primates)
%Create Newick-formatted character vector for benchmark
nwk15 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk15);
newNwkVect(16,1:nwkc)=nwk15;





%exp 3 dist seuclidean + upgma
S = std(descriptor,'omitnan');
distances=pdist(descriptor,'seuclidean',S);
UPGMAtree = seqlinkage(distances,'UPGMA',primates);
% h = plot(UPGMAtree,'orient','left');
% title('UPGMA Distance Tree of Primates');
% ylabel('Evolutionary distance');
%get canonical to identify the similarity
[Pointers, Distances, Names] = getcanonical(UPGMAtree);
%Create Newick-formatted character vector for benchmark
nwk16 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk16);
newNwkVect(17,1:nwkc)=nwk16;


UPGMAtree = seqlinkage(distances,'single',primates);
%Create Newick-formatted character vector for benchmark
nwk17 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk17);
newNwkVect(18,1:nwkc)=nwk17;

UPGMAtree = seqlinkage(distances,'complete',primates);
%Create Newick-formatted character vector for benchmark
nwk18 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk18);
newNwkVect(19,1:nwkc)=nwk18;

UPGMAtree = seqlinkage(distances,'weighted',primates);
%Create Newick-formatted character vector for benchmark
nwk19 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk19);
newNwkVect(20,1:nwkc)=nwk19;

UPGMAtree = seqlinkage(distances,'centroid',primates);
%Create Newick-formatted character vector for benchmark
nwk20 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk20);
newNwkVect(21,1:nwkc)=nwk20;

UPGMAtree = seqlinkage(distances,'median',primates);
%Create Newick-formatted character vector for benchmark
nwk21 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk21);
newNwkVect(22,1:nwkc)=nwk21;

% UPGMAtree1 = seqneighjoin(distances,'equivar',primates)
%Create Newick-formatted character vector for benchmark
nwk22 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk22);
newNwkVect(23,1:nwkc)=nwk22;

% UPGMAtree1 = seqneighjoin(distances,'firstorder',primates)
% %Create Newick-formatted character vector for benchmark
% nwk23 = getnewickstr(UPGMAtree1);
% [nwkr nwkc]=size(nwk23);
% newNwkVect(24,1:nwkc)=nwk23;




%exp 4 dist cityblock + upgma
distances=pdist(descriptor,'cityblock');
UPGMAtree = seqlinkage(distances,'UPGMA',primates);
%Create Newick-formatted character vector for benchmark
nwk24 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk24);
newNwkVect(25,1:nwkc)=nwk24;

UPGMAtree = seqlinkage(distances,'single',primates);
%Create Newick-formatted character vector for benchmark
nwk25 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk25);
newNwkVect(26,1:nwkc)=nwk25;

UPGMAtree = seqlinkage(distances,'complete',primates);
%Create Newick-formatted character vector for benchmark
nwk26 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk26);
newNwkVect(27,1:nwkc)=nwk26;

UPGMAtree = seqlinkage(distances,'weighted',primates);
%Create Newick-formatted character vector for benchmark
nwk27 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk27);
newNwkVect(28,1:nwkc)=nwk27;

UPGMAtree = seqlinkage(distances,'centroid',primates);
%Create Newick-formatted character vector for benchmark
nwk28 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk28);
newNwkVect(29,1:nwkc)=nwk28;

UPGMAtree = seqlinkage(distances,'median',primates);
%Create Newick-formatted character vector for benchmark
nwk29 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk29);
newNwkVect(30,1:nwkc)=nwk29;

UPGMAtree1 = seqneighjoin(distances,'equivar',primates)
%Create Newick-formatted character vector for benchmark
nwk30 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk30);
newNwkVect(31,1:nwkc)=nwk30;

UPGMAtree1 = seqneighjoin(distances,'firstorder',primates)
%Create Newick-formatted character vector for benchmark
nwk31 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk31);
newNwkVect(32,1:nwkc)=nwk31;



%exp 5 dist minkowski + upgma
distances=pdist(descriptor,'minkowski');
UPGMAtree = seqlinkage(distances,'UPGMA',primates);
%Create Newick-formatted character vector for benchmark
nwk32 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk32);
newNwkVect(33,1:nwkc)=nwk32;

UPGMAtree = seqlinkage(distances,'single',primates);
%Create Newick-formatted character vector for benchmark
nwk33 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk33);
newNwkVect(34,1:nwkc)=nwk33;

UPGMAtree = seqlinkage(distances,'complete',primates);
%Create Newick-formatted character vector for benchmark
nwk34 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk34);
newNwkVect(35,1:nwkc)=nwk34;

UPGMAtree = seqlinkage(distances,'weighted',primates);
%Create Newick-formatted character vector for benchmark
nwk35 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk35);
newNwkVect(36,1:nwkc)=nwk35;

UPGMAtree = seqlinkage(distances,'centroid',primates);
%Create Newick-formatted character vector for benchmark
nwk36 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk36);
newNwkVect(37,1:nwkc)=nwk36;

UPGMAtree = seqlinkage(distances,'median',primates);
%Create Newick-formatted character vector for benchmark
nwk37 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk37);
newNwkVect(38,1:nwkc)=nwk37;

UPGMAtree1 = seqneighjoin(distances,'equivar',primates)
%Create Newick-formatted character vector for benchmark
nwk38 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk38);
newNwkVect(39,1:nwkc)=nwk38;

UPGMAtree1 = seqneighjoin(distances,'firstorder',primates)
%Create Newick-formatted character vector for benchmark
nwk39 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk39);
newNwkVect(40,1:nwkc)=nwk39;





%exp 6 dist chebychev + upgma
distances=pdist(descriptor,'chebychev');
UPGMAtree = seqlinkage(distances,'UPGMA',primates);
%Create Newick-formatted character vector for benchmark
nwk40 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk40);
newNwkVect(41,1:nwkc)=nwk40;

UPGMAtree = seqlinkage(distances,'single',primates);
%Create Newick-formatted character vector for benchmark
nwk41 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk41);
newNwkVect(42,1:nwkc)=nwk41;

UPGMAtree = seqlinkage(distances,'complete',primates);
%Create Newick-formatted character vector for benchmark
nwk42 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk42);
newNwkVect(43,1:nwkc)=nwk42;

UPGMAtree = seqlinkage(distances,'weighted',primates);
%Create Newick-formatted character vector for benchmark
nwk43 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk43);
newNwkVect(44,1:nwkc)=nwk43;

UPGMAtree = seqlinkage(distances,'centroid',primates);
%Create Newick-formatted character vector for benchmark
nwk44 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk44);
newNwkVect(45,1:nwkc)=nwk44;

UPGMAtree = seqlinkage(distances,'median',primates);
%Create Newick-formatted character vector for benchmark
nwk45 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk45);
newNwkVect(46,1:nwkc)=nwk45;

UPGMAtree1 = seqneighjoin(distances,'equivar',primates)
%Create Newick-formatted character vector for benchmark
nwk46 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk46);
newNwkVect(47,1:nwkc)=nwk46;

UPGMAtree1 = seqneighjoin(distances,'firstorder',primates)
%Create Newick-formatted character vector for benchmark
nwk47 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk47);
newNwkVect(48,1:nwkc)=nwk47;






%exp 7 dist cosine + upgma
distances=pdist(descriptor,'cosine');
UPGMAtree = seqlinkage(distances,'UPGMA',primates);
%Create Newick-formatted character vector for benchmark
nwk48 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk48);
newNwkVect(49,1:nwkc)=nwk48;

UPGMAtree = seqlinkage(distances,'single',primates);
%Create Newick-formatted character vector for benchmark
nwk49 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk49);
newNwkVect(50,1:nwkc)=nwk49;

UPGMAtree = seqlinkage(distances,'complete',primates);
%Create Newick-formatted character vector for benchmark
nwk50 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk50);
newNwkVect(51,1:nwkc)=nwk50;

UPGMAtree = seqlinkage(distances,'weighted',primates);
%Create Newick-formatted character vector for benchmark
nwk51 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk51);
newNwkVect(52,1:nwkc)=nwk51;

UPGMAtree = seqlinkage(distances,'centroid',primates);
%Create Newick-formatted character vector for benchmark
nwk52 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk52);
newNwkVect(53,1:nwkc)=nwk52;

UPGMAtree = seqlinkage(distances,'median',primates);
%Create Newick-formatted character vector for benchmark
nwk53 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk53);
newNwkVect(54,1:nwkc)=nwk53;

UPGMAtree1 = seqneighjoin(distances,'equivar',primates)
%Create Newick-formatted character vector for benchmark
nwk54 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk54);
newNwkVect(55,1:nwkc)=nwk54;

UPGMAtree1 = seqneighjoin(distances,'firstorder',primates)
%Create Newick-formatted character vector for benchmark
nwk55 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk55);
newNwkVect(56,1:nwkc)=nwk55;





%exp 8 dist correlation + upgma
distances=pdist(descriptor,'correlation');
UPGMAtree = seqlinkage(distances,'UPGMA',primates);
%Create Newick-formatted character vector for benchmark
nwk56 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk56);
newNwkVect(57,1:nwkc)=nwk56;

UPGMAtree = seqlinkage(distances,'single',primates);
%Create Newick-formatted character vector for benchmark
nwk57 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk57);
newNwkVect(58,1:nwkc)=nwk57;

UPGMAtree = seqlinkage(distances,'complete',primates);
%Create Newick-formatted character vector for benchmark
nwk58 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk58);
newNwkVect(59,1:nwkc)=nwk58;

UPGMAtree = seqlinkage(distances,'weighted',primates);
%Create Newick-formatted character vector for benchmark
nwk59 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk59);
newNwkVect(60,1:nwkc)=nwk59;

UPGMAtree = seqlinkage(distances,'centroid',primates);
%Create Newick-formatted character vector for benchmark
nwk60 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk60);
newNwkVect(61,1:nwkc)=nwk60;

UPGMAtree = seqlinkage(distances,'median',primates);
%Create Newick-formatted character vector for benchmark
nwk61 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk61);
newNwkVect(62,1:nwkc)=nwk61;

UPGMAtree1 = seqneighjoin(distances,'equivar',primates)
%Create Newick-formatted character vector for benchmark
nwk62 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk62);
newNwkVect(63,1:nwkc)=nwk62;

UPGMAtree1 = seqneighjoin(distances,'firstorder',primates)
%Create Newick-formatted character vector for benchmark
nwk63 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk63);
newNwkVect(64,1:nwkc)=nwk63;





%exp 9 dist hamming + upgma
distances=pdist(descriptor,'hamming');
UPGMAtree = seqlinkage(distances,'UPGMA',primates);
%Create Newick-formatted character vector for benchmark
nwk64 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk64);
newNwkVect(65,1:nwkc)=nwk64;

UPGMAtree = seqlinkage(distances,'single',primates);
%Create Newick-formatted character vector for benchmark
nwk65 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk65);
newNwkVect(66,1:nwkc)=nwk65;

UPGMAtree = seqlinkage(distances,'complete',primates);
%Create Newick-formatted character vector for benchmark
nwk66 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk66);
newNwkVect(67,1:nwkc)=nwk66;

UPGMAtree = seqlinkage(distances,'weighted',primates);
%Create Newick-formatted character vector for benchmark
nwk67 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk67);
newNwkVect(68,1:nwkc)=nwk67;

UPGMAtree = seqlinkage(distances,'centroid',primates);
%Create Newick-formatted character vector for benchmark
nwk68 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk68);
newNwkVect(69,1:nwkc)=nwk68;

UPGMAtree = seqlinkage(distances,'median',primates);
%Create Newick-formatted character vector for benchmark
nwk69 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk69);
newNwkVect(70,1:nwkc)=nwk69;

UPGMAtree1 = seqneighjoin(distances,'equivar',primates)
%Create Newick-formatted character vector for benchmark
nwk70 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk70);
newNwkVect(71,1:nwkc)=nwk70;

UPGMAtree1 = seqneighjoin(distances,'firstorder',primates)
%Create Newick-formatted character vector for benchmark
nwk71 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk71);
newNwkVect(72,1:nwkc)=nwk71;





%exp 10 dist jaccard + upgma
distances=pdist(descriptor,'jaccard');
UPGMAtree = seqlinkage(distances,'UPGMA',primates);
%Create Newick-formatted character vector for benchmark
nwk72 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk72);
newNwkVect(73,1:nwkc)=nwk72;

UPGMAtree = seqlinkage(distances,'single',primates);
%Create Newick-formatted character vector for benchmark
nwk73 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk73);
newNwkVect(74,1:nwkc)=nwk73;

UPGMAtree = seqlinkage(distances,'complete',primates);
%Create Newick-formatted character vector for benchmark
nwk74 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk74);
newNwkVect(75,1:nwkc)=nwk74;

UPGMAtree = seqlinkage(distances,'weighted',primates);
%Create Newick-formatted character vector for benchmark
nwk75 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk75);
newNwkVect(76,1:nwkc)=nwk75;

UPGMAtree = seqlinkage(distances,'centroid',primates);
%Create Newick-formatted character vector for benchmark
nwk76 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk76);
newNwkVect(77,1:nwkc)=nwk76;

UPGMAtree = seqlinkage(distances,'median',primates);
%Create Newick-formatted character vector for benchmark
nwk77 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk77);
newNwkVect(78,1:nwkc)=nwk77;

UPGMAtree1 = seqneighjoin(distances,'equivar',primates)
%Create Newick-formatted character vector for benchmark
nwk78 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk78);
newNwkVect(79,1:nwkc)=nwk78;

UPGMAtree1 = seqneighjoin(distances,'firstorder',primates)
%Create Newick-formatted character vector for benchmark
nwk79 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk79);
newNwkVect(80,1:nwkc)=nwk79;




%exp 11 dist spearman + upgma
distances=pdist(descriptor,'spearman');
UPGMAtree = seqlinkage(distances,'UPGMA',primates);
%Create Newick-formatted character vector for benchmark
nwk80 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk80);
newNwkVect(81,1:nwkc)=nwk80;

UPGMAtree = seqlinkage(distances,'single',primates);
%Create Newick-formatted character vector for benchmark
nwk81 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk81);
newNwkVect(82,1:nwkc)=nwk81;

UPGMAtree = seqlinkage(distances,'complete',primates);
%Create Newick-formatted character vector for benchmark
nwk82 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk82);
newNwkVect(83,1:nwkc)=nwk82;

UPGMAtree = seqlinkage(distances,'weighted',primates);
%Create Newick-formatted character vector for benchmark
nwk83 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk83);
newNwkVect(84,1:nwkc)=nwk83;

UPGMAtree = seqlinkage(distances,'centroid',primates);
%Create Newick-formatted character vector for benchmark
nwk84 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk84);
newNwkVect(85,1:nwkc)=nwk84;

UPGMAtree = seqlinkage(distances,'median',primates);
%Create Newick-formatted character vector for benchmark
nwk85 = getnewickstr(UPGMAtree);
[nwkr nwkc]=size(nwk85);
newNwkVect(86,1:nwkc)=nwk85;

UPGMAtree1 = seqneighjoin(distances,'equivar',primates)
%Create Newick-formatted character vector for benchmark
nwk86 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk86);
newNwkVect(87,1:nwkc)=nwk86;

UPGMAtree1 = seqneighjoin(distances,'firstorder',primates)
%Create Newick-formatted character vector for benchmark
nwk87 = getnewickstr(UPGMAtree1);
[nwkr nwkc]=size(nwk87);
newNwkVect(88,1:nwkc)=nwk87;

%newick file write for reading in R file
fid = fopen('newWickTree.txt','w');
[nwr nwc]=size(newNwkVect);
for jj = 1 : nwr
  singleNwk=newNwkVect(jj,:);
  if jj~=24
      fprintf(fid, singleNwk);
      fprintf(fid, '\n');
  end
end
fclose(fid);

t=0;
t=t+1;