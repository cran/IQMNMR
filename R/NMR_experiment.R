NMR_experiment <-
function(
                         t.stimulated_metabolites_1,
                         s_data_x,
                         s_data_y,
                         variance_frequency,
                         variance_frequency_interval,
                         SNR,
                         TD,
                         SWH,
                         SFO1,
                         O1,
                         fid_save_path
                         ){

    signal_intensity_dss<-10*t.stimulated_metabolites_1[1,8]

    selected_metabolites<-t.stimulated_metabolites_1[which(t.stimulated_metabolites_1[,10]>0),]
    #rm(t.stimulated_metabolites_1, t.stimulated_metabolites_2)

    ID_signal_data<-selected_metabolites[-1,1] 

    #####  ============================== read space of gene  ================================
    #---------------------   1   read space of parameters-------------------------------------

    numbers_of_population<-2
    space_of_gene_x.ff<-s_data_x[which(match(s_data_x$ID,ID_signal_data)>0),]
    space_of_gene_y.ff<-s_data_y[ID_signal_data,]
    concentrations<-selected_metabolites[-1,10]
    sum_peak<-aggregate(space_of_gene_x.ff[,3], by=list(space_of_gene_x.ff[,1]), FUN=sum)
    total_nmr_signals<-((signal_intensity_dss/selected_metabolites[1,8])/selected_metabolites[1,10])*selected_metabolites[-1,8] #* (max_peak[,2]/sum_peak[,2])

    #------------------   2   select parameters randomly----------------------------------------------------------
    minimal_phase=-pi/2; maximal_phase=pi/2 
    minimal_damp_factor=5; maximal_damp_factor=20 
    chemical_shift_tep<-space_of_gene_x.ff[,2]
    amplitude_tep<-space_of_gene_x.ff[,3]
    damp_factor<-runif(length(chemical_shift_tep), min=minimal_damp_factor, max=maximal_damp_factor)
    phase<-runif(length(chemical_shift_tep), min=minimal_phase, max=maximal_phase)
    chemical_shift_y<-runif(length(space_of_gene_y.ff[,1]), min=-variance_frequency, max=variance_frequency)
    amplitude_y<-concentrations*total_nmr_signals   

    #------------------------------------  3   generate parameters----------------------------------------------------------
    # increasing the length of  population_chromosome_y
    y_times_chemical_shift.ff<-rep(chemical_shift_y,times=table(space_of_gene_x.ff[,1])[])
    y_times_amplitude.ff<-rep(amplitude_y,times=table(space_of_gene_x.ff[,1])[])
    sum_peak.ff<-rep(sum_peak[,2],times=table(space_of_gene_x.ff[,1])[])
    # fill in the parameters data ff
    #Chemical_Shifts
    chemical_shift1<-chemical_shift_tep*SFO1-abs(O1)+y_times_chemical_shift.ff
    chemical_shift<-chemical_shift1 + runif(length(chemical_shift1), min=-variance_frequency_interval, max=variance_frequency_interval)

    #amplitudes  
    amplitude<-(amplitude_tep/sum_peak.ff )*y_times_amplitude.ff 

    #----------------------------   4   generate output file-----------------------------------------------------
    stimulated_parameters<-matrix(0, nrow=(length(amplitude)+1), ncol=4)
    stimulated_parameters[-1,]<-cbind(chemical_shift, damp_factor, amplitude, phase)

    dss_parameters<-c(-abs(O1),20,signal_intensity_dss,0)
    stimulated_parameters[1,]<-dss_parameters
    rm(dss_parameters)
    min_signal_intensity<-stimulated_parameters[1,3]/stimulated_parameters[1,2]
    sd.noise1<-min_signal_intensity/SNR

    ### 04 FID
    chemical_shift_1<-c( stimulated_parameters[,1])
    amplitudes_1<-c(stimulated_parameters[,3])
    damping.factors_1<-c(stimulated_parameters[,2])
    phases_1<-c(stimulated_parameters[,4])
    rm(stimulated_parameters)

    fid1_function<-function(fk, dk, ak, pk, deltk, lengthk, sd.noise){
        t<-((1:lengthk)-1)*deltk
        para<-ak*exp(1i*pk)
        zk.para<-(-dk+1i*2*pi*fk)
        zk.Sli<-exp(outer(t, zk.para))
        zk.complex.mod<-ak*exp(1i*pk)
        fid.x<-zk.Sli%*%zk.complex.mod
        fid1<-Re(fid.x[])+rnorm(lengthk,mean=0,sd=sd.noise)
        fid2<-Im(fid.x[])+rnorm(lengthk,mean=0,sd=sd.noise)
        fid3<-t
        fid<-cbind(fid1, fid2, fid3)
        return(fid)
    }

    fid<-fid1_function(fk=chemical_shift_1, dk=damping.factors_1, ak=amplitudes_1, pk=phases_1, deltk=1/SWH, lengthk=TD/2, sd.noise=sd.noise1)
    pa<-cbind(chemical_shift_1, damping.factors_1, amplitudes_1, phases_1, amplitudes_1/damping.factors_1, round(amplitudes_1*1000))
    par<-pa[order(pa[,5]),]
    #write.table(par, file="par.csv", quote=TRUE, sep=";", row.names=FALSE)

    writefid<- function(fid.data.frame, writefile){
        wrfid<-c(1:(2*length(fid.data.frame[,1])))
        wrfid[seq(from=1,to=2*length(fid.data.frame[,1]) , by=2)]<-round(fid.data.frame[,1]*1000)
        wrfid[seq(from=2,to=2*length(fid.data.frame[,1]) , by=2)]<-round(fid.data.frame[,2]*1000)
        wrfid<-as.integer(wrfid)
        con.file<-file(writefile,"wb")
        writeBin(wrfid, writefile,  size=4,  endian="little")
        close(con.file)
    }

    writefid(fid.data.frame=fid, writefile=fid_save_path)
    
    mm1<-data.frame(max(fid[,3]),TD,SWH,SFO1,O1)
    names(mm1)<-c("AQ","TD","SWH","SFO1","O1")
    return(mm1)
}

