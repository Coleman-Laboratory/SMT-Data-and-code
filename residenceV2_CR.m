%% analysis the residence time

prompt = 'Please Define first cell in list that you want to Analyze:    ';
Start_CID = input(prompt); %set the first cell in the list

prompt = 'Please Define last cell in list that you want to Analyze:    ';
Last_CID = input(prompt); %set the last cell in the list

for cell_id=Start_CID:Last_CID
    ndwet=[];
% 
    acqu = sptana(cell_id).acqu;
    bleachrate= sptana(cell_id).bleachrate;
    ndwet = sptana(cell_id).ndwet;

    parg = 2;
    binsize = 2*acqu;
    [pop_1, pop_2 , tau_1, tau_2]=pdf_fitCR(ndwet, binsize, sptana(cell_id).bleachrate, parg, 1);
    
    pop1(cell_id,1) = pop_1;
    pop2(cell_id,1) = pop_2;
   
    tau1(cell_id,1) = tau_1;
    tau2(cell_id,1) = tau_2;
     
end


  