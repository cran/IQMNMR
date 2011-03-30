select_metabolites <-
function(n_metabolites,
                             s.lists_metabolites
                             ){
    number_H=9
    contentration_DSS=5
    maximal_concentration=30
    minimal_concentration=0.5
    ID_signal_data<-sort(sample(s.lists_metabolites$ID,n_metabolites))
    signal_data<-s.lists_metabolites
    
    concentration<-c(rep(0,dim(s.lists_metabolites)[1]))
    concentration2<-runif(n=n_metabolites, min=minimal_concentration, max=maximal_concentration)
    concentration[c(ID_signal_data)]<-concentration2
       
    signal_data<-as.data.frame(cbind(signal_data,concentration))
    signal_data_dss<-as.data.frame(list(0,"DSS", signal_data[1,3], signal_data[1,4], signal_data[1,5], signal_data[1,6], "zg", number_H, "C6 H16 O3 S Si",  contentration_DSS))
    names(signal_data_dss)<-names(signal_data)
    signal_data<-rbind(signal_data_dss,signal_data)
    selected_metabolites<-signal_data
    return(selected_metabolites)
}

