NMR_spectra <-
function(SWH,
                      AQ,
                      TD,
                      fidfile,
                      SFO1,
                      O1
                      ){
    #read fid of bruker
    readfid<-function(r.SWH, r.pulse.time, r.number_fidbin, r.fidfile){
         require(ff)
         tpulse.time<-r.pulse.time   
         to.read = file(r.fidfile,"rb")  
         signal<-readBin(to.read, ,what="int",size=4, n=r.number_fidbin, endian = "little")
         close(to.read)
         td <- length(signal)
         rawR <- signal[seq(from = 1, to = td, by = 2)]
         rawI <- signal[seq(from = 2, to = td, by = 2)]
         mediar<-mean(as.integer(rawR[c((3*length(rawR)/4):length(rawR))]),na.rm = TRUE)
         mediai<--mean(as.integer(rawI[c((3*length(rawR)/4):length(rawR))]),na.rm = TRUE)
         rawR<-rawR-mediar
         rawI<-rawI-mediai
         fid.raw<-ff(1:((td/2)*2),vmode="double",dim=c((td/2),2), update=TRUE, overwrite=TRUE)
         fid.raw[,1]<-c(rawR)  
         fid.raw[,2]<-c(rawI)  
         rm("rawR","rawI")
         return(fid.raw)
    }


    fftshow<-function(f.swh, re.fid, im.fid, location.TMS, f.SFO1){
        require(ff)

        N<-length(re.fid); dwell.time<-1/f.swh; f.max<-f.swh
        dwell.f<-f.max/(4*N)  

        f.set.positive<-c(c(0:(2*N-1))*dwell.f)
        f.set.negative<-rev(-c(c(1:(2*N-1))*dwell.f))
        f.set<-c(f.set.negative,f.set.positive)

        rm("f.set.positive","f.set.negative")

        fid<-c(re.fid)+c(im.fid)*1i

        zero.re.fid<-c(c(fid),c(rep(0,3*N)))

        rm("fid")

        fft.shift<-function (x){
            lon<-length(x)
            y <- numeric(lon-1)
            y[c(1:(lon/2-1))] <- x[(lon/2 + 1):(lon-1)]
            y[c((lon/2):(lon-1))] <- x[1:(lon/2)]
            return(y)
        }
    
        re.freq<-fft.shift(fft(zero.re.fid))

        m<-length(re.freq)

        rm("zero.re.fid","fft.shift")
 
        freq.date<-ff(c(1:(m*4)), vmode="double",dim=c(m,4), update=TRUE, overwrite=TRUE)
        freq.date[,1]<-Re(re.freq)
        freq.date[,2]<-Im(re.freq)
        freq.date[,3]<-f.set

        freq.date[,4]<- ((f.set-location.TMS)/SFO1)

        return(freq.date)
    }

    fid1<-readfid(r.SWH=SWH, r.pulse.time=AQ, r.number_fidbin=TD, r.fidfile=fidfile)
    

    freq1<-fftshow(f.swh=SWH, re.fid=fid1[,1], im.fid=fid1[,2], location.TMS=-abs(O1), f.SFO1=SFO1)
    
    
    par(mfcol=c(2,1))
    plot(fid1[,1],type="l", ylab="re",xlab="TD")
    plot(freq1[,1]~freq1[,4],type="l", xlim=c(max(freq1[,4]),min(freq1[,4])),ylim=c(min(freq1[,1]),max(freq1[,1])),ylab="re",xlab="ppm")
    
    result<-list(fid1,freq1)
    names(result)<-c("spectrum_of_time_domain","spectrum_of_frequency_domain")
    return(result)
}

