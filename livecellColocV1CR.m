%%determine colocalization points from two color data collection
% input file is the sptana file where two color data must be entered into spt_process_v15CR as a listing of cells in
% the first color channel followed by a list of the corresponded cells in
% the second color. For example, an experiment where 10 cells were
% collected in two colors must have the first ten rows of sptana as cells
% 1-10 in color 1 followed by rows 11-20 of cells in color 2
% The program will find nearest molecules in channel/color 1 that are within a user specified minimum colocalization distance from molecules in channel/color 2 in pixels.
% This program saves the following values, percentages of molecules
% colocalized in each cell (perColoc), percentages of molecules temporally
% colocalized in each cell (perClcT),aggregated arrival/departure times of channel/color 2 molecules with respect to channel/color 1 (ADBag), aggregated colocalization times (ClcTag),
% and a cell structure with with a table containing all molecular IDs (from non-segmented binding events), arrival departure frames, XY position, arrival/departure times of channel/color 2 with respect to channel/color 1 (FinalMoldata))

clear all

ADBag=[];
ClcTag=[];
FinalMoldata=[];
MolBfil=[];

[FileName, Pathname] = uigetfile({'.mat'},'Please grab sptana file');
file=[Pathname, FileName];
sptana = importdata(file);
%%
prompt = 'Please input total number of cells with two color data from sptana file:    '; %input number of cells
nm = input(prompt); % allows program to grab the second color data from corresponding cell

prompt = 'Please Define the acquisition time in seconds:   ';
acqu = input(prompt); % set the acquistion time + mechanical delay in seconds, used 1.34 for PCNA/PolD two color experiment

prompt = 'Please enter maximum colocalization distance in pixels:    '; % number of pixels to be considered as co-localized
coloc_min= input(prompt); 

prompt = 'Filter out any molecules less than this time (seconds), enter 0 for no filter:    '; % Allows user to filter out transient interactions, used 1.35 for PCNA and PolD experiments
filter= input(prompt); 


for n=1:nm;
    MlocA=[];
    MlocB=[];
    
     MlocA(:,1)=sptana(n).nucl_summary(:,8);
     MlocA(:,2)=sptana(n).nucl_summary(:,9);
     MlocB(:,1)=sptana(n+nm).nucl_summary(:,8);
     MlocB(:,2)=sptana(n+nm).nucl_summary(:,9);
     MolAdwet=sptana(n).ndwet;
     MolBdwet=sptana(n+nm).ndwet;
    FidxAdwet=find(MolAdwet>=filter);
    FidxBdwet=find(MolBdwet>=filter);
    MColoc=[];
    coloc=[];
    colocd=[];
    MolAnz=[];
    MolBnz=[];
    MolA=MlocA(FidxAdwet,1:2);
    MolB=MlocB(FidxBdwet,1:2);
    MolBfil(n,1)=length(MolB);
    MolAnz=MolA;
    molBco=[];
    molAco=[];
    Idxco=[];
    MolBnz=MolB;
    IDX = knnsearch(MolAnz(:,1:2),MolBnz,'K',4); %gather 4 closest nearest neighbor
    MColoc=MolAnz(IDX(:,1),:); %defines which nearest neighbor you are analyzing
    coloc=MolBnz - MColoc;
    colocd=((coloc(:,1).^2+coloc(:,2).^2).^0.5); %calculates the distance between molecules in pixels
    colocidx=colocd < coloc_min;
    molBco=MolB(colocidx,1:2);
    IDXco=IDX(colocidx);
    molAco=MolA(IDXco,1:2);
    perColoc(n,1)=length(molAco)/length(MolB); %calculates percent of molecules throughout move that are colocalized
    rtA=[];
    rtB=[];
    colocTA=[];
    colocTB=[];
    rtA=knnsearch(MlocA(:,1:2),molAco,'K',1);
    rtB=knnsearch(MlocB(:,1:2),molBco,'K',1);
    colocTA=sptana(n).nucl_summary(rtA,:);
    colocTB=sptana(n+nm).nucl_summary(rtB,:);
    StartB=[];
    EndB=[];
    ClcT=[];
    MolIDA=[];
    MolIDB=[];
    ResiClcTA=[];
    ResiClcTB=[];
    for tx=1:length(colocTB(:,1));
        testIndB=[];
        testIndA=[];
        testIndTC=[];
        testIndB(1:1000,1)=zeros;
        testIndB(colocTB(tx,2):colocTB(tx,3),1)=ones;
        testIndA(1:1000,1)=zeros;
        testIndA(colocTA(tx,2):colocTA(tx,3),1)=ones;
        testIndTC=testIndB.*testIndA;
        
        if length(find(testIndTC(:,1)==1))>=1;
            if colocTB(tx,2)>=2;
                StartB(tx,1)=(colocTB(tx,2)-colocTA(tx,2))*acqu;
                EndB(tx,1)=(colocTB(tx,3)-colocTA(tx,3))*acqu;
                ClcT(tx,1)=(length(find(testIndTC(:,1)==1)))*acqu;
                ResiClcTA(tx,1)=(colocTA(tx,3)-colocTA(tx,2)+1)*acqu;
                ResiClcTB(tx,1)=(colocTB(tx,3)-colocTB(tx,2)+1)*acqu;
                MolIDA(tx,1:3)=colocTA(tx,1:3);
                MolIDB(tx,1:3)=colocTB(tx,1:3);
                MolIDA(tx,4:5)=colocTA(tx,8:9);
                MolIDB(tx,4:5)=colocTB(tx,8:9);
                
            end
        end
        

    end
    FinalMol=[];
    FinalMol=[MolIDA MolIDB StartB EndB ClcT ResiClcTA ResiClcTB];
    stBFx=[];
    ClcTFx=[];
    StartBF=[];
    EndBF=[];
    MolIDAFx=[];
    MolIDBFx=[];
    MolIDAF=[];
    MolIDBF=[];
    
    if length(ClcT)>=1;
        ClcTFx=find(ClcT(:,1)==0);
        ClcTF=ClcT;
        ClcTF(ClcTFx)=[];
        StartBF=StartB;
        StartBF(ClcTFx)=[];
        EndBF=EndB;
        EndBF(ClcTFx)=[];
        MolIDAF=MolIDA;
        MolIDAF(ClcTFx,:)=[];
        MolIDBF=MolIDB;
        MolIDBF(ClcTFx,:)=[];
        FinalMol(ClcTFx,:)=[];
        
    end
    FinalMoldata{n,1}=FinalMol; %creates Final table for each cell with mol id, arrival/departure time, xy position of colocalized molecules, colocalized time, residence time molecule in Channel/color 1, residence time molecule in Channel/color 2 
    perClcT(n,1)=length(StartBF)/length(MolB); %calculates percent of colocalizing molecules that temporally overlap
    ADBtemp=[];
    ClcTtemp=[];
    if length(EndBF)>=1
        ADBtemp=[];
        ADBtemp(:,1)=StartBF(:,1);
        ADBtemp(:,2)=EndBF(:,1);
        ADB{n,1}=ADBtemp;
        ClcTtemp=[];
        ClcTtemp=ClcTF(:,1);
    end
    ADBag=cat(1,ADBag,ADBtemp); %generates an aggregated list of arrival (column 1) and departure (column 2) channel/color 2 molecules with respect to arrival/departure of molecules in Channel/color 1
    ClcTag=cat(1,ClcTag,ClcTtemp); %generates an aggregated list of colocalization times for molecules within all cells analyzed
end
%% save the processed data in file

uisave({'sptana','perClcT','ADBag','ClcTag','FinalMoldata','perColoc'});