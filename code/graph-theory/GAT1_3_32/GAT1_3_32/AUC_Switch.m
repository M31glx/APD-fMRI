AUCorFDA = input('Type 1 for AUC or 2 for FDA analysis: ');

switch  AUCorFDA
    
    case 1 %AUC analysis
        
        AUC_type = input('Type 1 for global, 2 for regional, or 3 for hub measures: ');

        switch AUC_type

            case 1%global

                AUC_Analysis_NetMes

            case 2%regional

                AUC_Analysis_NetReg

            case 3%hubs

                AUC_NetHubs      
        end
        
    case 2%FDA
        
        FDA_type = input('Type 1 for global, 2 for regional, or 3 for hub measures: ');

        switch FDA_type

            case 1%global

                FDA_Analysis_NetMes

            case 2%regional

                FDA_Analysis_NetReg

            case 3%hubs

                FDA_NetHubs      
        end
end

        
