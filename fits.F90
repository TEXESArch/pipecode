
      subroutine capword(word)
!       capitalize an 8-character word

        character(8) word

        do i = 1,8
          ichari = ichar(word(i:i))
          if((ichari>96).and.(ichari<=122)) then
            ichari = ichari-32
            word(i:i) = char(ichari)
          endif
        enddo

        return
      end


      subroutine capline(line)
!       capitalize the first 8 characters of a line

        character(80) line

        do i = 1,8
          ichari = ichar(line(i:i))
          if((ichari>96).and.(ichari<=122)) then
            ichari = ichari-32
            line(i:i) = char(ichari)
          endif
        enddo

        return
      end


      subroutine fithcomm(keyword,comment,iunit)
!	write a comment or history fits header line to iunit

        character(8) keyword
        character(60) comment
        character(80) line

        call capword(keyword)
        lcom = lastchar(comment)
        lcom = min(lcom,60)
        if(lcom>0) then
          line = keyword//"  "//comment(1:lcom)
        endif
        write(iunit,'(a80)') line

        return
      end


      subroutine fithchar(keyword,string,comment,iunit)
!	write a fits header line for a character type keyword

        character(8) keyword
        character(60) string,comment
        character(80) line

        call capword(keyword)
        last = lastchar(string)
        line = keyword//"= '"//string(1:last)//"'"
        last = max((last+12),30)
        lcom = lastchar(comment)
        lcom = min(lcom,(76-last))
        if(lcom>0) then
          line = line(1:last)//"  / "//comment(1:lcom)
        endif
        write(iunit,'(a80)') line

        return
      end


      function lastchar(string)
!	finds the last non-blank character in a string

        character(60) string

        ic = 61
  100   ic = ic-1
        if((ic>0).and.(string(ic:ic)==' ')) goto 100
        lastchar = ic

        return

      end


      subroutine fithlog(keyword,lval,comment,iunit)
!	write a fits header line for a logical type keyword

        character(8) keyword
        character(60) comment
        character(80) line
        logical lval

        call capword(keyword)
        write(unit=line,fmt='(a8,"= ",l20)') keyword,lval
        last = 30
        lcom = lastchar(comment)
        lcom = min(lcom,(76-last))
        if(lcom>0) then
          line = line(1:last)//"  / "//comment(1:lcom)
        endif
        write(iunit,'(a80)') line

        return
      end


      subroutine fithint(keyword,ival,comment,iunit)
!	write a fits header line for a integer type keyword

        character(8) keyword
        character(60) comment
        character(80) line

        call capword(keyword)
        write(unit=line,fmt='(a8,"= ",i20)') keyword,ival
        last = 30
        lcom = lastchar(comment)
        lcom = min(lcom,(76-last))
        if(lcom>0) then
          line = line(1:last)//"  / "//comment(1:lcom)
        endif
        write(iunit,'(a80)') line

        return
      end


      subroutine fithreal(keyword,val,comment,iunit)
!	write a fits header line for a real type keyword

        character(8) keyword
        character(60) comment
        character(80) line

        call capword(keyword)
        write(unit=line,fmt='(a8,"= ",1p,e20.6)') keyword,val
        last = 30
        lcom = lastchar(comment)
        lcom = min(lcom,(76-last))
        if(lcom>0) then
          line = line(1:last)//"  / "//comment(1:lcom)
        endif
        write(iunit,'(a80)') line

        return
      end



      subroutine fitshd(oldfile,ierr)
!       read oldfile.hd, convert to fits format, and write to fits.hd

        use ius
        use modes
        use consts

        character(32) oldfile,oldhd
        character(80) line,val
        character(60) object,feature,obsmode,instmode,objtype, &
                      piname,pid,weather,warning
        character(60) comment,NOC
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical leftopen
        real nodpa,lores,kmirror,krot

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
                piname,pid,weather,warning
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl

        NOC = ' '
        ra = 15.*12+34/4.+56.78/240.
        equinox = 2000.0
        epoch = 2000.0
        pltscl = 0.
        dvpix = 0.

!       open old and fits header files

        print '(" Enter name of file for header copy")'
        read '(a24)', oldfile
        namelen = index(oldfile,'.hd')-1
        if(namelen<=0) namelen = index(oldfile,' ')-1
        oldhd = oldfile(:namelen)//'.hd'
        open(unit=iurawh,file=oldhd)

        inquire(unit=iufith,opened=leftopen)
        if(leftopen) close(unit=iufith)
        open(unit=iufith,file='fits.hd')

!       read, write, and parse lines

        ierr = 0
        gain = 1.0
        obstime = 1.0
        jorder = 0
  100   read(unit=iurawh,fmt='(a79)',err=190,end=200) line
        comment = NOC
        if(line(1:8)=='object  ') then
          read(unit=line,fmt='(10x,a16)',err=190,end=190) object
          call fithchar('OBJECT  ',object,comment,iufith)
        elseif(line(1:8)=='feature ') then
          read(unit=line,fmt='(10x,a16)',err=190,end=190) feature
          comment = 'spectral feature'
          call fithchar('FEATURE ',feature,comment,iufith)
        elseif(line(1:8)=='obsmode ') then
          read(unit=line,fmt='(10x,a16)',err=190,end=190) obsmode
          if(index(obsmode,'chop-nod')>0) then
            comment = 'OBJECT'
          elseif(index(obsmode,'nod')>0) then
            comment = 'OBJECT'
          elseif(index(obsmode,'chop')>0) then
            comment = 'OBJECT'
          elseif(index(obsmode,'stare')>0) then
            comment = 'ENG'
          elseif(index(obsmode,'dark')>0) then
            comment = 'ENG'
          elseif(index(obsmode,'flat')>0) then
            comment = 'FLAT'
          elseif(index(obsmode,'map')>0) then
            comment = 'OBJECT'
          elseif(index(obsmode,'scan')>0) then
            comment = 'OBJECT'
          elseif(index(obsmode,'fowler')>0) then
            comment = 'OBJECT'
          else
            comment = 'ENG'
          endif
          call fithchar('OBSMODE ',obsmode,comment,iufith)
        elseif(line(1:8)=='instmode') then
          read(unit=line,fmt='(10x,a16)',err=190,end=190) instmode
          if(index(instmode,'hi-med')>0) then
            comment = 'echelon x echelle'
          elseif(index(instmode,'med')>0) then
            comment = 'echelle only'
          elseif(index(instmode,'hi-lo')>0) then
            comment = 'echelon x lo-res'
          elseif(index(instmode,'lo')>0) then
            comment = 'lo-res only'
          elseif(index(instmode,'cam')>0) then
            comment = 'lo-res face-on'
          elseif(index(instmode,'pup')>0) then
            comment = 'pupil viewing camera'
          endif
          call fithchar('INSTMODE',instmode,comment,iufith)
        elseif(line(1:8)=='waveno0 ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) waveno0
          if(waveno0<=350.) waveno0 = 1000.
          comment = 'planned central wavenumber'
          call fithreal('WAVENO0 ',waveno0,comment,iufith)
        elseif(line(1:8)=='order   ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) horder
          jorder = nint(horder)
          comment = 'cross-dispersion grating order number'
          call fithint('ORDER   ',jorder,comment,iufith)
        elseif(line(1:8)=='temp    ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) temper
          comment = 'ambient temperature'
          if((temper>(-15.)).and.(temper<25.)) then
            temp = temper+273.16
          elseif((temper<250.).or.(temper>300.)) then
            print '(" Implausible temperature:",1p,1e10.2)', temper
            temp = 273.16
            comment = 'substituted for implausible value'
          else
            temp = temper
          endif
          call fithreal('TEMPER  ',temp,comment,iufith)
        elseif(line(1:8)=='nod     ') then
          if(index(obsmode,'scan')<=0) then
            iche = index(line,'E')
            ichw = max(iche,index(line,'W'))
            ichn = index(line,'N')
            if((iche<12).or.(ichn.lt.iche+2)) goto 190
            val = line(10:iche-1)
            read(unit=val,fmt=*,err=190,end=190) dra
            val = line(ichw+1:ichn-1)
            read(unit=val,fmt=*,err=190,end=190) ddec
            if(ddec==0.) then
              dist = dra
            elseif(dra==0.) then
              dist = ddec
            else
              dist = sqrt(dra**2+ddec**2)
            endif
            call fithreal('NOD     ',dist,line(11:70),iufith)
          endif
        elseif(line(1:8)=='telstep ') then
          if(index(obsmode,'scan')>0) then
            iche = index(line,'E')
            ichn = index(line,'N')
            if((iche<12).or.(ichn.lt.iche+2)) goto 190
            val = line(10:iche-1)
            read(unit=val,fmt=*,err=190,end=190) dra
            val = line(iche+1:ichn-1)
            read(unit=val,fmt=*,err=190,end=190) ddec
            if(ddec==0.) then
              dist = dra
            elseif(dra==0.) then
              dist = ddec
            else
              dist = sqrt(dra**2+ddec**2)
            endif
            call fithreal('TELSTEP ',dist,line(11:70),iufith)
          endif
        elseif((line(1:8)=='inttime ') &
                .or.(line(1:8)=='frametim')) then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) frtime
          if(frtime<=0.) frtime = 1.0
          comment = 'integration time / frame'
          call fithreal('FRAMETIM',frtime,comment,iufith)
        elseif(line(1:8)=='obstime ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) obstime
          if(obstime<=0.) obstime = 1.0
          comment = 'calculated observation time'
          call fithreal('OBSTIME ',obstime,comment,iufith)
        elseif(line(1:8)=='tottime ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) tottime
          if(tottime<=0.) tottime = 1.0
          comment = 'calculated total clock time'
          call fithreal('TOTTIME ',tottime,comment,iufith)
        elseif(line(1:8)=='echelle ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) echelle
          echelle = echelle-xdmr0
          if((modeinst==MMED).or.(modeinst.eq.MHIMED)) &
            xdr = tan(echelle/DEGRAD)
          comment = 'echelle grating angle (deg)'
          call fithreal('ECHELLE ',echelle,comment,iufith)
        elseif(line(1:8)=='lores   ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) lores
          lores = lores-xdlr0
          if((modeinst==MLOW).or.(modeinst.eq.MHILOW)) &
            xdr = tan(lores/DEGRAD)
          comment = 'lo-res grating angle (deg)'
          call fithreal('LORES   ',lores,comment,iufith)
        elseif(line(1:8)=='kmirror ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) kmirror
          comment = 'K-mirror angle'
          call fithreal('KMIRROR ',kmirror,comment,iufith)
        elseif(line(1:8)=='slit    ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) slit
          comment = 'slit wheel angle'
          call fithreal('SLIT    ',slit,comment,iufith)
        elseif(line(1:8)=='filter  ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) filter
          comment = 'filter wheel angle'
          call fithreal('FILTER  ',filter,comment,iufith)
        elseif(line(1:8)=='gain    ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) gain
          if(gain==0.) gain = 1.0
          comment = 'electronic gain'
          call fithreal('GAIN    ',gain,comment,iufith)
        elseif(line(1:8)=='pixelwd ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) pixelwd
!       may have to convert pixelwd to cm here
          comment = 'pixel width (um?)'
          call fithreal('PIXELWD ',pixelwd,comment,iufith)
        elseif(line(1:8)=='dvpix   ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) dvpix
          comment = 'Doppler shift per pixel (km/s)'
          call fithreal('DVPIX   ',dvpix,comment,iufith)
        elseif(line(1:8)=='pltscl  ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) pltscl
          if(abs(slitpa-180.)<0.1) pltscl = -abs(pltscl)
          comment = 'platescale (arcsec/pixel)'
          call fithreal('PLTSCL  ',pltscl,comment,iufith)
        elseif(line(1:8)=='nframe  ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nframe
          comment = 'number of frames coadded in hardware'
          call fithint('NFRAME  ',nframe,comment,iufith)
        elseif(line(1:8)=='nsum    ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nsum
          comment = 'number of frames coadded in software'
          call fithint('NSUM    ',nsum,comment,iufith)
        elseif(line(1:8)=='nwrite  ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nwrite
          comment = 'number of frames written per nod phase'
          call fithint('NWRITE  ',nwrite,comment,iufith)
        elseif(line(1:8)=='nnod    ') then
          if((modeobs==MSCAN).or.(modeobs.eq.MMAP)) then
            read(unit=line,fmt='(10x,i8)',err=190,end=190) nscan
            comment = 'number of scans'
            call fithint('NSCAN   ',nscan,comment,iufith)
          else
            read(unit=line,fmt='(10x,i8)',err=190,end=190) nnod
            comment = 'number of nod pairs'
            call fithint('NNOD    ',nnod,comment,iufith)
          endif
        elseif((line(1:8)=='nsteps  ').or.(line(1:8).eq.'npoints ')) &
          then
          if((modeobs==MSCAN).or.(modeobs.eq.MMAP)) then
            if(line(1:8)=='nsteps   ') then
              read(unit=line,fmt='(10x,i8)',err=190,end=190) nnod
            else
              read(unit=line,fmt='(10x,f10.0)',err=190,end=190) points
              nnod = points
            endif
          endif
          comment = 'number of points in scan'
          call fithint('NPOINTS ',nnod,comment,iufith)
        elseif(line(1:8)=='nsky    ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nsky
          comment = 'number of extra points on sky'
          call fithint('NSKY    ',nsky,comment,iufith)
        elseif(line(1:8)=='nspec   ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nspec
          ny = nspec
          comment = 'array size in spectral dimentsion'
          call fithint('NSPEC   ',nspec,comment,iufith)
        elseif(line(1:8)=='nspat   ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nspat
          comment = 'array size in spatial dimentsion'
          call fithint('NSPAT   ',nspat,comment,iufith)
        elseif(line(1:8)=='bitpix  ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nbitpix
        elseif(line(1:8)=='time    ') then
          read(unit=line,fmt='(10x,i2,1x,i2,1x,f3.0)', &
            err=191,end=191) ih,im,s
          ut = ih+im/60.+s/3600.
          comment = 'UTC ?'
          call fithchar('TIME-OBS',line(11:70),comment,iufith)
        elseif(line(1:8)=='date-obs') then
          read(unit=line,fmt='(10x,i2,1x,i2,1x,i2)', &
            err=109,end=191) iyr,imo,idy
          goto 110
  109     read(unit=line,fmt='(10x,i4,1x,i2,1x,i2)', &
            err=191,end=191) iyr,imo,idy
          iyr = iyr-2000
  110     date = iyr+(imo+idy/100.)/100.
          call fithchar('DATE-OBS',line(11:70),NOC,iufith)
          if(imo<=2) then
            jyr = iyr+2001
            jmo = imo+6
            jdy = idy
          elseif(imo<=6) then
            jyr = iyr+2001
            jmo = imo+6
            jdy = min(idy,30)
          elseif(imo==8) then
            jyr = iyr+2002
            jmo = imo-6
            jdy = min(idy,28)
          else
            jyr = iyr+2002
            jmo = imo-6
            jdy = min(idy,30)
          endif
          write(unit=comment,fmt='(i4,"-",i2.2,"-",i2.2)') jyr,jmo,jdy
          call fithchar('DATE-REL',comment,NOC,iufith)
        elseif(line(1:8)=='ra      ') then
          read(unit=line,fmt='(10x,i2,1x,i2,1x,f5.2)', &
            err=191,end=191) ih,im,s
          ra = 15.*ih+im/4.+s/240.
          call fithreal('RA      ',ra,line(11:21),iufith)
        elseif(line(1:8)=='dec     ') then
          if(line(11:11)=='-') then
            read(unit=line,fmt='(10x,i3,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            dec = id-im/60.-s/3600.
          elseif(line(11:11)=='+') then
            read(unit=line,fmt='(11x,i2,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            dec = id+im/60.+s/3600.
          else
            read(unit=line,fmt='(10x,i2,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            dec = id+im/60.+s/3600.
          endif
          call fithreal('DEC     ',dec,line(11:21),iufith)
        elseif(line(1:8)=='azimuth ') then
          if(line(11:11)=='-') then
            read(unit=line,fmt='(10x,i4,1x,i2,1x,f5.2)', &
              err=192,end=191) id,im,s
            az = id-im/60.-s/3600.
          elseif(line(11:11)=='+') then
            read(unit=line,fmt='(11x,i3,1x,i2,1x,f5.2)', &
              err=192,end=191) id,im,s
            az = id+im/60.+s/3600.
          else
            read(unit=line,fmt='(10x,i3,1x,i2,1x,f5.2)', &
              err=192,end=191) id,im,s
            az = id+im/60.+s/3600.
          endif
          call fithreal('AZIMUTH ',az,line(11:21),iufith)
        elseif(line(1:8)=='elevatio') then
          if(line(11:11)=='-') then
            read(unit=line,fmt='(10x,i3,1x,i2,1x,f5.2)', &
              err=193,end=191) id,im,s
            el = id-im/60.-s/3600.
          elseif(line(11:11)=='+') then
            read(unit=line,fmt='(11x,i2,1x,i2,1x,f5.2)', &
              err=193,end=191) id,im,s
            el = id+im/60.+s/3600.
          else
            read(unit=line,fmt='(10x,i2,1x,i2,1x,f5.2)', &
              err=193,end=191) id,im,s
            el = id+im/60.+s/3600.
          endif
          call fithreal('ELEVATIO',el,line(11:21),iufith)
        elseif(line(1:8)=='equinox ') then
          read(unit=line,fmt='(10x,f7.2)',err=191,end=191) equinox
          if(equinox<1900.) then
            read(unit=line,fmt='(10x,i4)',err=191,end=191) ieq
            equinox = ieq
          endif
          epoch = equinox
          call fithreal('EQUINOX ',equinox,NOC,iufith)
        elseif(line(1:8)=='ha      ') then
          read(unit=line,fmt='(10x,f7.4)',err=194,end=191) ha
          call fithreal('HA      ',ha,NOC,iufith)
        elseif(line(1:8)=='airmass ') then
          read(unit=line,fmt='(10x,f7.3)',err=191,end=191) airmass
          call fithreal('AIRMASS ',airmass,NOC,iufith)
        elseif(line(1:8)=='instrpa ') then
          read(unit=line,fmt='(10x,f7.3)',err=191,end=191) slitpa
          call fithreal('INSTRPA ',slitpa,NOC,iufith)
        elseif(line(1:8)=='instraa ') then
          read(unit=line,fmt='(10x,f7.3)',err=191,end=191) aa
          call fithreal('INSTRAA ',aa,NOC,iufith)
        elseif(line(1:8)=='telescop') then
          if(index(line,'McD')>0) then
            efl0 = 12.*270.
          elseif(index(line,'IRTF')>0) then
            efl0 = 12.*300.
          elseif(index(line,'Gemini')>0) then
            if(date<(6.4)) then
              efl0 = 10.6*800.
            else
              efl0 = 10.7*800.
            endif
          else
            print '(" Unknown telescope")'
            efl0 = 12.*300.
          endif
          call fithchar('TELESCOP',line(11:70),NOC,iufith)
        elseif(line(1:7)=='instrum') then
          comment = 'TEXES'
          call fithchar('INSTRUME',comment,NOC,iufith)
        elseif((line(1:7)=='nx     ') &
          .or.(line(1:7)=='ny     ') &
          .or.(line(1:7)=='nz     ')) then
          continue
        else
          call fithchar(line(1:8),line(11:70),NOC,iufith)
        endif
        goto 100
  192   continue
          if(line(11:11)=='-') then
            read(unit=line,fmt='(10x,i3,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            az = id-im/60.-s/3600.
          elseif(line(11:11)=='+') then
            read(unit=line,fmt='(11x,i2,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            az = id+im/60.+s/3600.
          else
            read(unit=line,fmt='(10x,i2,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            az = id+im/60.+s/3600.
          endif
          call fithreal('AZIMUTH ',az,line(11:21),iufith)
        goto 100
  193   continue
          if(line(11:11)=='-') then
            read(unit=line,fmt='(10x,i2,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            el = id-im/60.-s/3600.
          elseif(line(11:11)=='+') then
            read(unit=line,fmt='(11x,i1,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            el = id+im/60.+s/3600.
          else
            read(unit=line,fmt='(10x,i1,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            el = id+im/60.+s/3600.
          endif
          call fithreal('ELEVATIO',el,line(11:21),iufith)
  194   continue
          if(line(11:11)=='-') then
            read(unit=line,fmt='(10x,i3,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            ha = id-im/60.-s/3600.
          elseif(line(11:11)=='+') then
            read(unit=line,fmt='(11x,i2,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            ha = id+im/60.+s/3600.
          else
            read(unit=line,fmt='(10x,i2,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            ha = id+im/60.+s/3600.
          endif
          call fithreal('HA      ',ha,line(11:21),iufith)
        goto 100

  190   ierr = ierr+1
  191   print '(" Error reading line from ",a16/a80)', oldhd,line
        goto 100

  200   if(pltscl==0.) then
          print '(" Enter platescale: ",$)'
          read *, pltscl
          comment = 'platescale (arcsec/pixel)'
          call fithreal('PLTSCL  ',pltscl,comment,iufith)
        endif
        if(dvpix==0.) then
          print '(" Enter dv/pixel: ",$)'
          read *, dvpix
          dvpix = -abs(dvpix)
          comment = 'Doppler shift per pixel (km/s)'
          call fithreal('DVPIX   ',dvpix,comment,iufith)
        endif

        return

      end



      subroutine fitsrdhd(iufits,irec,ierr)
!       read a fits header from an opened fits file

        use dims
        use modes
        use consts

        character(80) line
        character(30) val
        character(2880) fitsrec
        character(60) object,feature,obsmode,instmode,objtype, &
                      piname,pid,weather,warning
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        real nodpa,lores,kmirror,krot
        character(8) yn
        logical ifobj

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
                piname,pid,weather,warning
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /tort0/ hrfl0,xdfl0,xdmrdgr,xdlrdgr,xdmr0,xdlr0,fred0

        ierr = 0
        gain = 1.0
        obstime = 1.0
        jorder = 0
        ra = 15.*12+34/4.+56.78/240.
        equinox = 2000.0
        naxis = 1
        naxis1 = 1
        naxis2 = 1
        naxis3 = 1

!       read a 2880 byte record from fitsfile
        irec = 0
  100   irec = irec+1
        read(unit=iufits,rec=irec,err=190) fitsrec

!	interpret record as 36 80-character lines
        do ichar = 1,2880,80
          line = fitsrec(ichar:ichar+79)
          if(line(11:11)==char(39)) then
            last = index(line(12:30),char(39))+10
            if(last>11) then
              val = line(12:last)
            else
              val = line(12:41)
            endif
          else
            val = line(11:40)
          endif
!	first 10 parameters not currently stored in common
          if(line(1:8)=='NAXIS   ') then
            read(unit=val,fmt=*,err=190,end=190) naxis
          elseif(line(1:8)=='NAXIS1  ') then
            read(unit=val,fmt=*,err=190,end=190) naxis1
          elseif(line(1:8)=='NAXIS2  ') then
            read(unit=val,fmt=*,err=190,end=190) naxis2
          elseif(line(1:8)=='NAXIS3  ') then
            read(unit=val,fmt=*,err=190,end=190) naxis3
          elseif(line(1:8)=='CRVAL2  ') then
            read(unit=val,fmt=*,err=190,end=190) crval2
          elseif(line(1:8)=='CRPIX2  ') then
            read(unit=val,fmt=*,err=190,end=190) crpix2
          elseif(line(1:8)=='CDELT2  ') then
            read(unit=val,fmt=*,err=190,end=190) cdelt2
          elseif(line(1:8)=='CTYPE2  ') then
            if(index(val,'VELO')>0) then
              dvpix = cdelt2
              pix0 = crpix2
              vel0 = crval2
            endif
          elseif(line(1:8)=='CRVAL3  ') then
            read(unit=val,fmt=*,err=190,end=190) crval3
          elseif(line(1:8)=='CRPIX3  ') then
            read(unit=val,fmt=*,err=190,end=190) crpix3
          elseif(line(1:8)=='CDELT3  ') then
            read(unit=val,fmt=*,err=190,end=190) cdelt3
          elseif(line(1:8)=='CTYPE3  ') then
            if(index(val,'VELO')>0) then
              dvpix = cdelt3
              pix0 = crpix3
              vel0 = crval3
            endif
          elseif(line(1:8)=='OBJECT  ') then
            read(unit=val,fmt='(a16)',err=190,end=190) object
          elseif(line(1:8)=='FEATURE ') then
            read(unit=val,fmt='(a16)',err=190,end=190) feature
          elseif(line(1:8)=='OBSMODE ') then
            read(unit=val,fmt='(a16)',err=190,end=190) obsmode
            if(index(obsmode,'chop-nod')>0) then
              modeobs = MCHOPNOD
            elseif(index(obsmode,'nod')>0) then
              modeobs = MNOD
            elseif(index(obsmode,'chop')>0) then
              modeobs = MCHOP
            elseif(index(obsmode,'stare')>0) then
              modeobs = MSTARE
            elseif(index(obsmode,'dark')>0) then
              modeobs = MSTARE
            elseif(index(obsmode,'flat')>0) then
              modeobs = MFLAT
            elseif(index(obsmode,'map')>0) then
              modeobs = MMAP
            elseif(index(obsmode,'scan')>0) then
              modeobs = MSCAN
            elseif(index(obsmode,'cell')>0) then
              modeobs = MCELL
            elseif(index(obsmode,'fowler')>0) then
              modeobs = MNOD
            else
!              modeobs = MUNKNOWN
              print '("*** Unknown obs mode ",a20)', obsmode
            endif
            ifobj = (modeobs/=MFLAT)
          elseif(line(1:8)=='INSTMODE') then
            read(unit=val,fmt='(a16)',err=190,end=190) instmode
            if(index(instmode,'hi-med')>0) then
              modeinst = MHIMED
            elseif(index(instmode,'med')>0) then
              modeinst = MMED
            elseif(index(instmode,'hi-lo')>0) then
              modeinst = MHILOW
            elseif(index(instmode,'lo')>0) then
              modeinst = MLOW
            elseif(index(instmode,'cam')>0) then
              modeinst = MCAMERA
            elseif(index(instmode,'pup')>0) then
              modeinst = MCAMERA
            else
!              modeinst = MUNKNOWN
              print '("*** Unknown inst mode ",a20)', instmode
            endif
            if((modeinst==MMED).or.(modeinst.eq.MHIMED)) then
              xddgr = xdmrdgr
              xdr = tan(echelle/DEGRAD)
            elseif((modeinst==MLOW).or.(modeinst.eq.MHILOW)) then
              xddgr = xdlrdgr
              xdr = tan(lores/DEGRAD)
            endif
          elseif(line(1:8)=='WAVENO0 ') then
            read(unit=val,fmt=*,err=190,end=190) waveno0
            if(waveno0<=400.) waveno0 = 1000.
          elseif(line(1:8)=='WNO0    ') then
            read(unit=val,fmt=*,err=190,end=190) wno0
            if(wno0<=400.) wno0 = 1000.
          elseif(line(1:8)=='ORDER   ') then
            read(unit=val,fmt=*,err=190,end=190) horder
            jorder = nint(horder)
          elseif(line(1:8)=='TEMPER  ') then
            read(unit=val,fmt=*,err=190,end=190) temper
            if((temper>(-15.)).and.(temper<25.)) then
              temp = temper+273.16
            elseif((temper<250.).or.(temper>300.)) then
              print '(" Implausible temperature:",1p,1e10.2)', temper
              temp = 273.16
            else
              temp = temper
            endif
          elseif(line(1:8)=='NOD     ') then
            read(unit=val,fmt=*,err=190,end=190) dist
          elseif(line(1:8)=='TELSTEP ') then
            read(unit=val,fmt=*,err=190,end=190) dist
          elseif(line(1:8)=='FRAMETIM') then
            read(unit=val,fmt=*,err=190,end=190) frtime
            if(frtime<=0.) frtime = 1.0
          elseif(line(1:8)=='OBSTIME ') then
            read(unit=val,fmt=*,err=190,end=190) obstime
            if(obstime<=0.) obstime = 1.0
          elseif(line(1:8)=='TOTTIME ') then
            read(unit=val,fmt=*,err=190,end=190) tottime
            if(tottime<=0.) tottime = 1.0
          elseif(line(1:8)=='ECHELLE ') then
            read(unit=val,fmt=*,err=190,end=190) echelle
            echelle = echelle-xdmr0
            if((modeinst==MMED).or.(modeinst.eq.MHIMED)) &
              xdr = tan(echelle/DEGRAD)
          elseif(line(1:8)=='LORES   ') then
            read(unit=val,fmt=*,err=190,end=190) lores
            lores = lores-xdlr0
            if((modeinst==MLOW).or.(modeinst.eq.MHILOW)) &
              xdr = tan(lores/DEGRAD)
          elseif(line(1:8)=='KMIRROR ') then
            read(unit=val,fmt=*,err=190,end=190) kmirror
          elseif(line(1:8)=='SLIT    ') then
            read(unit=val,fmt=*,err=190,end=190) slit
          elseif(line(1:8)=='FILTER  ') then
            read(unit=val,fmt=*,err=190,end=190) filter
          elseif(line(1:8)=='GAIN    ') then
            read(unit=val,fmt=*,err=190,end=190) gain
            if(gain==0.) gain = 1.0
          elseif(line(1:8)=='PIXELWD ') then
            read(unit=val,fmt=*,err=190,end=190) pixelwd
!       may have to convert pixelwd to cm here
          elseif(line(1:8)=='DVPIX   ') then
            read(unit=val,fmt=*,err=190,end=190) dvpix
          elseif(line(1:8)=='PLTSCL  ') then
            read(unit=val,fmt=*,err=190,end=190) pltscl
            if(abs(slitpa-180.)<0.1) pltscl = -abs(pltscl)
          elseif(line(1:8)=='NFRAME  ') then
            read(unit=val,fmt=*,err=190,end=190) nframe
          elseif(line(1:8)=='NSUM    ') then
            read(unit=val,fmt=*,err=190,end=190) nsum
          elseif(line(1:8)=='NWRITE  ') then
            read(unit=val,fmt=*,err=190,end=190) nwrite
          elseif(line(1:8)=='NNOD    ') then
            read(unit=val,fmt=*,err=190,end=190) nnod
            if(nnod>mp) then
              print '(" nnod =",i3," > mp =",i3, &
                "  Reading only through mp")', nnod,mp
              nnod = mp
            endif
          elseif(line(1:8)=='NSCAN   ') then
            read(unit=val,fmt=*,err=190,end=190) nscan
          elseif((line(1:8)=='NSTEPS  ').or.(line(1:8).eq.'NPOINTS ')) &
            then
            if((modeobs==MSCAN).or.(modeobs.eq.MMAP)) then
              if(line(1:8)=='nsteps   ') then
                read(unit=val,fmt=*,err=190,end=190) nnod
              else
                read(unit=val,fmt=*,err=190,end=190) points
                nnod = points
              endif
            endif
          elseif(line(1:8)=='NSKY    ') then
            read(unit=val,fmt=*,err=190,end=190) nsky
          elseif(line(1:8)=='NSPEC   ') then
            read(unit=val,fmt=*,err=190,end=190) nspec
            ny = nspec
          elseif(line(1:8)=='NSPAT   ') then
            read(unit=val,fmt=*,err=190,end=190) nspat
          elseif(line(1:8)=='BITPIX  ') then
            read(unit=val,fmt=*,err=190,end=190) nbitpix
          elseif(line(1:4)=='DATE') then
            read(unit=val,fmt=*,iostat=ierr) date
            if(ierr/=0) then
              read(unit=val,fmt='(i4,1x,i2,1x,i2)',iostat=ierr) &
                iyr,imo,idy
              iyr = iyr-2000
              if(ierr/=0) then
                ierr = 0
                read(unit=val,fmt='(i2,1x,i2,1x,i2)',err=190,end=190) &
                  iyr,imo,idy
              endif
              ierr = 0
              date = iyr+(imo+idy/100.)/100.
            endif
            if(date>2000.) date = date-2000.
            if((date>30.).or.(date<0.)) then
              print '(" date bad ",f12.4,a30)', date,val
              stop
            endif
          elseif(line(1:8)=='RA      ') then
            read(unit=val,fmt=*,err=190,end=190) ra
          elseif(line(1:8)=='DEC     ') then
            read(unit=val,fmt=*,err=190,end=190) dec
          elseif(line(1:8)=='AZIMUTH ') then
            read(unit=val,fmt=*,err=190,end=190) az
          elseif(line(1:8)=='ELEVATIO') then
            read(unit=val,fmt=*,err=190,end=190) el
          elseif(line(1:8)=='HA      ') then
            read(unit=val,fmt=*,err=190,end=190) ha
          elseif(line(1:8)=='AIRMASS ') then
            read(unit=val,fmt=*,err=191,end=191) airmass
          elseif(line(1:8)=='INSTRPA ') then
            read(unit=val,fmt=*,err=191,end=191) slitpa
          elseif(line(1:8)=='TELESCOP') then
            if(index(line,'McD')>0) then
              efl0 = 12.*270.
            elseif(index(line,'IRTF')>0) then
              efl0 = 12.*300.
            elseif(index(line,'Gemini')>0) then
              if(date<(6.4)) then
                efl0 = 10.6*800.
              else
                efl0 = 10.7*800.
              endif
            else
              print '(" Unknown telescope")'
              efl0 = 12.*300.
            endif
          elseif(line(1:8)=='INSTRUM ') then
            continue
          elseif(line(1:8)=='END     ') then
            goto 120
          endif
        enddo
        goto 100

  120   xnkbr = sqrt(2.3613-3.115e+04/waveno0**2-5.86e+08/waveno0**4)
!       index at 10um
        xnkbr0 = 1.5263
        if(date<(1.03)) then
!       ZnSe focal reducer
          fred = 1.93*((waveno0-320.)/(waveno0-300.))/(680./700.)
          nc = 3
        else
!       KBr focal reducer
          pow = .525*(xnkbr-1.)/(xnkbr0-1.)
          fred = 1./(1.-pow)
          nc = nnod
        endif
        if(irtf) then
          flfore = 47.25/(xnkbr-1.)
          flfor0 = 47.25/(xnkbr0-1.)
          flcdte = 240.
          d01 = 125.
          d12 = 350.
!         xmfore = (dfore-flfore)/flfore
!         xmfor0 = (dfore-flfor0)/flfor0
          xmfore = flfore*flcdte/((d01-flfore)*(d12-flcdte)-d01*flfore)
          xmfor0 = flfor0*flcdte/((d01-flfor0)*(d12-flcdte)-d01*flfor0)
          forerat = xmfore/xmfor0
        else
          forerat = 1.
        endif
!       adjust fred to get fls right
        fred = 1.01*fred
!       not sure when this changed
        if(date>(16.0)) fred = 1.01*fred
        hrfl = hrfl0/fred
        xdfl = xdfl0/fred
!       efl = efl0*forerat/fred
        efl = efl0/(forerat*fred)
        pltscl = pixelwd/(efl*4.848e-06)
        slitwid = slitval(slit,date)/(efl0*4.848e-06)
        omegap = pltscl*slitwid*(4.848e-06)**2
        if((modeinst==MHIMED).or.(modeinst.eq.MMED)) then
          iorder = nint(2.0*xddgr*sin(echelle/DEGRAD)*waveno0)
          if(jorder*(iorder-jorder)/=0) &
            print '(" Warning: header and calculated order disagree:", &
              2i4/3es10.3)', jorder,iorder,xddgr,echelle,waveno0
          wno = iorder/(2.0*xddgr*sin(echelle/DEGRAD))
          if(abs(wno-waveno0)>(waveno0/500.)) then
            print '(" Warning: echelle and waveno0 disagree", &
              2f10.2)', echelle,waveno0
            print '(" Derived order, waveno: ",i2,f10.2)', iorder,wno
            write(iupl,'(" Warning: echelle and waveno0 disagree")')
            if((abs(abs(wno0)-wno)<(wno/100.)) &
              .or.(abs(abs(wno0)-waveno0)<(waveno0/100.))) then
              print '(" Setting waveno0 = wno0")'
              waveno0 = abs(wno0)
              if(abs(wno-waveno0)>(waveno0/500.)) then
                print '(" Adjust mr0? (y)",$)'
                read '(a8)', yn
                if(index(yn,'n')==0) then
                  xdorder = iorder
                  xdang0 = DEGRAD*asin(xdorder/(2.0*xddgr*waveno0))
                  xdmr0 = echelle+xdmr0-xdang0
                  echelle = xdang0
                  print '(" mr0 = ",f10.3)', xdmr0
                  write(iupl,'(" mr0 = ",f10.3)') xdmr0
                endif
              endif
            elseif(abs(waveno0-wno)<(waveno0/100.)) then
              print '(" Setting waveno0 = calculated")'
              waveno0 = wno
            else
              print '("^G*** xdr may be wrong in tort")'
            endif
          endif
        elseif((modeinst==MHILOW).or.(modeinst.eq.MLOW)) then
          iorder = nint(2.0*xddgr*sin(lores/DEGRAD)*waveno0)
          if(jorder*(iorder-jorder)/=0) &
            print '(" Warning: header and calculated order disagree:", &
              2i4/3es10.3)', jorder,iorder,xddgr,lores,waveno0
          wno = iorder/(2.0*xddgr*sin(lores/DEGRAD))
          if(abs(wno-waveno0)>(waveno0/200.)) then
            print '(" Warning: lores and waveno0 disagree", &
              3f10.2)', lores,wno,waveno0
            write(iupl,'(" Warning: lores and waveno0 disagree")')
            print '(" Derived order, waveno: ",i2,f10.2)', iorder,wno
            if((abs(abs(wno0)-wno)<(wno/50.)) &
              .or.(abs(abs(wno0)-waveno0)<(waveno0/50.))) then
              waveno0 = abs(wno0)
              print '(" Setting waveno0 = wno0")'
              if(abs(wno-abs(wno0))>(wno/100.)) then
                print '(" Adjust lr0? (y)",$)'
                read '(a8)', yn
                if(index(yn,'n')==0) then
                  xdorder = iorder
                  xdang0 = DEGRAD*asin(xdorder/(2.0*xddgr*waveno0))
                  xdlr0 = lores+xdlr0-xdang0
                  lores = xdang0
                  print '(" lr0 = ",f10.3)', xdlr0
                  write(iupl,'(" lr0 = ",f10.3)') xdlr0
                endif
              endif
            elseif(abs(waveno0-wno)<(waveno0/50.)) then
              waveno0 = wno
              print '(" Setting waveno0 = wno")'
            else
              print '("^G*** xdr may be wrong in tort")'
            endif
          endif
        endif
!       I think waveno0 is more reliable than echelle for calculating xdr
!       wno = iorder/(2.0*xddgr*sin(echelle/DEGRAD))
        if(modeinst/=MCAMERA) then
          sinang = iorder/(2.0*xddgr*waveno0)
          xdr = sinang/sqrt(1.-sinang**2)
        endif

        if(index(obsmode,'fowler')>0) then
          nframe = nframe/2
          frtime = (52./77.)*frtime/nframe
        endif
        beamtime = frtime*nframe*nsum*nwrite

        if(ifobj.and.(modeobs==MSCAN)) then
          if(date<(1.0)) then
            nsky = 0
          elseif(date<(1.07)) then
            nsky = (nnod+9)/10
            nsky = max0(3,nsky)
          endif
          nnod = nnod+nsky
        endif
        if(nnod>mp) then
          print '(" nnod =",i3," > mp =",i3, &
                "  Arrays will overflow")', nnod,mp
          write(iupl,'(" nnod =",i3," > mp =",i3, &
                "  Arrays will overflow")') nnod,mp
          ierr = -1
          return
        endif

        if(naxis==3) then
          if(naxis3<256) then
            nu = naxis1
            nv = naxis2+4
            nw = naxis3
          else
            nu = naxis3
            nv = naxis1+4
            nw = naxis2
          endif
        else
          ns = naxis1
          nt = naxis2+4
        endif
        if(wno0==1000.) wno0 = waveno0/(1.0+vel0/2.998e+05)
        return
  
  190   ierr = ierr+1
  191   print '(" Error reading fits header line ",2i6/a80)', &
          irec,ichar,line
        if(ierr>7) return
        goto 100

        return

      end



      subroutine fitsread(iufits,irec,txfile,ierr)
!       read fits data and write out to txfile

        use dims
        use modes
        use consts

        character(32) txfile
        character(60) object,feature,obsmode,instmode,objtype, &
                      piname,pid,weather,warning
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        real nodpa,lores,kmirror,krot
        logical doflip,doflop

        common /byteflip/doflip,doflop
        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
                piname,pid,weather,warning
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        dimension idata(mdata),idatx(mx)
        dimension fdata(mdata),fdatx(mx)

        open(unit=iurawd,file=txfile,access='direct',recl=4*nx,iostat=ierr)
        if(ierr/=0) then
          print '(" Error opening ",a32)', txfile
          stop
        endif
          

        if((modeobs==MSTARE).or.(modeobs==MCELL)) then
          nread = nwrite*nnod
        elseif(modeobs==MFLAT) then
          nread = nwrite*nnod
        elseif((modeobs==MNOD).or.(modobs==MCHOP)) then
          nread = 2*nwrite*nnod
        elseif(modeobs==MCHOPNOD) then
          nread = 4*nwrite*nnod
        elseif(modeobs==MSCAN) then
          nread = nwrite*nnod*nscan
        elseif(modeobs==MMAP) then
          nread = 2*nwrite*nnod*nscan
        endif
        ndat = mdata
        idat = ndat
        jrec = 0
        do iread = 1,nread
          do iy = 1,ny
            jrec = jrec+1
            do ix = 1,nx
              idat = idat+1
              if(nbitpix>0) then
                if(idat>ndat) then
                  irec = irec+1
                  read(unit=iufits,rec=irec,err=190) idata
                  if(doflip) call flipi4(idata,ndat)
                  idat = 1
                endif
                idatx(ix) = idata(idat)
              else
                if(idat>ndat) then
                  irec = irec+1
                  read(unit=iufits,rec=irec,err=190) fdata
                  if(doflip) call flipr4(fdata,fdata,ndat)
                  idat = 1
                endif
                fdatx(ix) = fdata(idat)
              endif
            enddo
            if(nbitpix>0) then
              if(doflip) call flipi4(idatx,nx)
              write(unit=iurawd,rec=jrec) idatx
            else
              if(doflip) call flipr4(fdatx,fdatx,nx)
              write(unit=iurawd,rec=jrec) fdatx
            endif
          enddo
          lread = iread
        enddo
        close(unit=iurawd)
        return

  190   print '(" Error reading fits data",i4," of",i4" frames read")', &
          lread,nread
        ierr = -1
        close(unit=iurawd)

        return
      end



      subroutine rawstart(fitsfile,nz,ierr)
!	open a raw fits file and write the primary header from fits.hd
!	leave open to write data

        use ius
        use modes

        character(48) fitsfile
        character(80) line
        logical baddata,notdone
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /nn/ nx,ny,nc,norder,ns,nt
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf

        ierr = 0
        if(.not.dofits) return
        call fitsopen(iufraw,fitsfile,ierr)
        if(ierr/=0) then
          print '(" Error opening ",a48)', fitsfile
          return
        endif

        naxis = 3
        write(line,'("SIMPLE  = ",l20)') .true.
        call fitsline(iufraw,line)
        write(line,'("BITPIX  = ",i20)') 32
        call fitsline(iufraw,line)
        write(line,'("NAXIS   = ",i20)') naxis
        call fitsline(iufraw,line)
        write(line,'("NAXIS1  = ",i20)') nx
        call fitsline(iufraw,line)
        write(line,'("NAXIS2  = ",i20)') ny
        call fitsline(iufraw,line)
        write(line,'("NAXIS3  = ",i20)') nz
        call fitsline(iufraw,line)
        write(line,'("EXTEND  = ",l20)') .true.
        call fitsline(iufraw,line)
        write(line,'("COMMENT   ", &
          "Raw data is stored as nx x ny x nz =",i4,"x",i4,"x",i4, &
          " integer array")') nx,ny,nz
        call fitsline(iufraw,line)
        close(unit=iufith)
        open(unit=iufith,file='fits.hd')
        notdone = .true.
        do
          read(iufith,'(a80)',iostat=ierr) line
          if(ierr/=0) exit
          if(index(line,'PIPELINE')>0) notdone = .false.
          if(notdone) call fitsline(iufraw,line)
        enddo
        backspace(unit=iufith)

        write(line,'("ORIGIN  = ''University of Texas''")')
        call fitsline(iufraw,line)
        write(line,'("END")')
        call fitsline(iufraw,line)
        call fitslpad(iufraw)
        ierr = 0

        return
      end


      subroutine fithedit(iufits,line)
!	replace a line in a fits file header

        character(80) line
        character(2880) fitsrec
        logical ifopen

        inquire(unit=iufits,opened=ifopen)
        if(.not.ifopen) then
          print '("*** Attempt to edit unopened fits file ",i4)', ifits
          return
        endif

        irec = 1
        read(unit=iufits,rec=irec) fitsrec
        ichline = index(fitsrec,line(1:10))
        if(ichline<1) then
          print '(1x,a10," not found in fits header")', line(1:10)
        else
          fitsrec(ichline:ichline+79) = line
          write(unit=iufits,rec=irec) fitsrec
          print '(" Editing fits header line")'
          print '(a80)', line
        endif

        return
      end


      subroutine rawxtend(nx,ny,iz)
!	write the extension header for a raw frame

        use ius
        use modes

        character(80) line
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf

        if(.not.dofits) return

        naxis = 2

        line = "XTENSION= 'IMAGE   '"
        call fitsline(iufraw,line)
        write(line,'("BITPIX  = ",i20)') 32
        call fitsline(iufraw,line)
        write(line,'("NAXIS   = ",i20)') naxis
        call fitsline(iufraw,line)
        write(line,'("NAXIS1  = ",i20)') nx
        call fitsline(iufraw,line)
        write(line,'("NAXIS2  = ",i20)') ny
        call fitsline(iufraw,line)
        write(line,'("PCOUNT  = ",i20)') 0
        call fitsline(iufraw,line)
        write(line,'("GCOUNT  = ",i20)') 1
        call fitsline(iufraw,line)
        line = "EXTNAME = 'RAWFRAME'"
        call fitsline(iufraw,line)
        if(iz>0) then
          write(line,'("EXTVER  = ",i20)') iz
          call fitsline(iufraw,line)
        endif

        write(line,'("END")')
        call fitsline(iufraw,line)
        call fitslpad(iufraw)

        return
      end


      subroutine redstart(nx,ny,nz,fitsfile,iunit,ierr)
!	open a reduced fits file and write the primary header
!	copying header lines from iunit = iufith or iufish

        use ius
        use modes

        character(48) fitsfile
        character(80) line
        logical baddata,hdopen
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf

        ierr = 0
        if(.not.dofits) return
        call fitsopen(iufred,fitsfile,ierr)
        if(ierr/=0) print '("*** Error opening ",a48)', fitsfile

        if(nx<=0) then
          naxis = 0
        elseif(nz>1) then
          naxis = 3
        elseif(ny>1) then
          naxis = 2
        elseif(nx>1) then
          naxis = 1
        else
          naxis = 0
        endif

        write(line,'("SIMPLE  = ",l20)') .true.
        call fitsline(iufred,line)
        write(line,'("BITPIX  = ",i20)') -32
        call fitsline(iufred,line)
        write(line,'("NAXIS   = ",i20)') naxis
        call fitsline(iufred,line)
        if(naxis>0) then
          write(line,'("NAXIS1  = ",i20)') nx
          call fitsline(iufred,line)
          if(naxis>1) then
            write(line,'("NAXIS2  = ",i20)') ny
            call fitsline(iufred,line)
            if(naxis>2) then
              write(line,'("NAXIS3  = ",i20)') nz
              call fitsline(iufred,line)
            endif
          endif
        endif
        write(line,'("EXTEND  = ",l20)') .true.
        call fitsline(iufred,line)
        if((modeobs==MMAP).or.(modeobs.eq.MSCAN)) then
          if(naxis/=3) then
            print '("*** redstart was called with bad dimensions")'
            print '(" nx,ny,nz,naxis = ",4i8)', nx,ny,nz,naxis
          endif
          write(line,'("COMMENT   ", &
            "The primary data unit contains a spectral-spatial-spatial")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            " data cube")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
           "The primary HDU is followed by a table extension")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
           "(EXTNAME = EXTRACTED)")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "giving the observed vacuum wavenumbers,")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "flux through the slit, estimated flux uncertainty,")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "estimated atmo transmission = (black-sky)/black")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "and the observed vacuum wavelengths")')
          call fitsline(iufred,line)
!          write(line,'("COMMENT   ", &
!           "The primary HDU is followed by four 1-D IMAGE extensions")')
!          call fitsline(iufred,line)
!          write(line,'("COMMENT   ", &
!            "giving the wavenumber scale (EXTNAME = WAVENUMBER),")')
!          call fitsline(iufred,line)
!          write(line,'("COMMENT   ", &
!            "the extracted flux spectrum (may be zero) (FLUX),")')
!          call fitsline(iufred,line)
!          write(line,'("COMMENT   ", &
!            "the calculated 1-sigma noise in the flux (NOISE),")')
!          call fitsline(iufred,line)
!          write(line,'("COMMENT   ", &
!            "and the telluric transmission spectrum calculated from")')
!          call fitsline(iufred,line)
!          write(line,'("COMMENT   ", &
!            "(black-sky)/black (EXTNAME = ATMO)")')
!          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "After the table extension are two 2-D extensions")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "giving the sky frame which has been subtracted from the")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "scan frames (EXTNAME = SKY),")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "and the calculated noise frame (EXTNAME = SCAN-NOISE)")')
          call fitsline(iufred,line)
          write(line,'("CTYPE1  = ''spectral pixel''")')
          call fitsline(iufred,line)
          write(line,'("CTYPE2  = ''pixel along slit''")')
          call fitsline(iufred,line)
          write(line,'("CTYPE3  = ''pixel along scan''")')
          call fitsline(iufred,line)
        elseif(modeinst<MCAMERA) then
          if(naxis/=2) &
            print '("*** redstart was called with bad dimensions")'
          write(line,'("COMMENT   ", &
            "The primary HDU contains a 2-D spectral-spatial image")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
           "The primary HDU is followed by a table extension")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
           "(EXTNAME = EXTRACTED)")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "giving the observed vacuum wavenumbers,")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "flux through the slit, estimated flux uncertainty,")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "estimated atmo transmission = (black-sky)/black")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "and the observed vacuum wavelengths")')
          call fitsline(iufred,line)
!          write(line,'("COMMENT   ", &
!           "The primary HDU is followed by four 1-D IMAGE extensions")')
!          call fitsline(iufred,line)
!          write(line,'("COMMENT   ", &
!            "giving the wavenumber scale (EXTNAME = WAVENUMBER),")')
!          call fitsline(iufred,line)
!          write(line,'("COMMENT   ", &
!            "the extracted flux spectrum (EXTNAME = FLUX),")')
!          call fitsline(iufred,line)
!          write(line,'("COMMENT   ", &
!            "the calculated 1-sigma noise in the flux (NOISE),")')
!          call fitsline(iufred,line)
!          write(line,'("COMMENT   ", &
!            "and the telluric transmission spectrum calculated from")')
!          call fitsline(iufred,line)
!          write(line,'("COMMENT   ", &
!            "(black-sky)/black (EXTNAME = ATMO)")')
!          call fitsline(iufred,line)
          write(line,'("CTYPE1  = ''spectral pixel''")')
          call fitsline(iufred,line)
          write(line,'("CTYPE2  = ''pixel along slit''")')
          call fitsline(iufred,line)
        else
          if(naxis/=2) &
            print '("*** redstart was called with bad dimensions")'
          write(line,'("COMMENT   ", &
            "The primary HDU contains a 2-D spatial-spatial image")')
          call fitsline(iufred,line)
          write(line,'("COMMENT   ", &
            "It may be followed by a 2-D spatial-spatial noise frame")')
          call fitsline(iufred,line)
        endif

!	copy iunit to iufred
        close(unit=iunit)
        if(iunit==iufith) then
          open(unit=iunit,file='fits.hd')
        elseif(iunit==iufish) then
          open(unit=iunit,file='fitsum.hd')
        else
          print '("*** Invalid temporary fits header unit",i4)', iunit
          ierr = 8
          return
        endif
        do
          read(iunit,'(a80)',iostat=ierr) line
          if(ierr/=0) exit
          call fitsline(iufred,line)
        enddo
        if(iunit==iufith) then
          close(unit=iunit)
        else
          backspace(unit=iunit)
        endif
        ierr = 0

        write(line,'("ORIGIN  = ''University of Texas''")')
        call fitsline(iufred,line)
        write(line,'("END")')
        call fitsline(iufred,line)
        call fitslpad(iufred)

        return
      end


      subroutine redxtend(iufitf,nx,ny,iz,itype)
!	write the header for an extension to a reduced fits file

        use ius
        character(80) line
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf

        if(.not.dofits) return
        if(ny>1) then
          naxis = 2
        elseif(nx>1) then
          naxis = 1
        else
          naxis = 0
        endif

        line = "XTENSION= 'IMAGE   '"
        call fitsline(iufred,line)
        write(line,'("BITPIX  = ",i20)') -32
        call fitsline(iufred,line)
        write(line,'("NAXIS   = ",i20)') naxis
        call fitsline(iufred,line)
        write(line,'("NAXIS1  = ",i20)') nx
        call fitsline(iufred,line)
        if(naxis>1) then
          write(line,'("NAXIS2  = ",i20)') ny
          call fitsline(iufred,line)
        endif
        write(line,'("PCOUNT  = ",i20)') 0
        call fitsline(iufred,line)
        write(line,'("GCOUNT  = ",i20)') 1
        call fitsline(iufred,line)
        write(line,'("CTYPE1  = ''spectral pixel''")')
        call fitsline(iufred,line)
        if(itype==1) then
          line = "EXTNAME = 'WAVENUMBER'" &
            //"  / wavenumbers for corresponding data values"
          call fitsline(iufred,line)
          line = "BUNIT   = 'cm^-1'"
          call fitsline(iufred,line)
        elseif(itype==2) then
          line = "EXTNAME = 'FLUX'" &
            //"  / fluxes extracted from 2-D data"
          call fitsline(iufred,line)
          line = "BUNIT   = 'Jy  / (may be inaccurate)'"
          call fitsline(iufred,line)
        elseif(itype==3) then
          line = "EXTNAME = 'NOISE'" &
            //"  / calculated 1-sigma noise in extracted fluxes"
          call fitsline(iufred,line)
          line = "BUNIT   = 'Jy  / (may be inaccurate)'"
          call fitsline(iufred,line)
        elseif(itype==4) then
          line = "EXTNAME = 'ATMO'" &
            //"  / atmospheric transmission calculated from black-sky"
          call fitsline(iufred,line)
          line = "BUNIT   = 'none'"
          call fitsline(iufred,line)
        elseif(itype==(-1)) then
          line = "EXTNAME = 'SKY'" &
            //"  / sky frame which was subtracted from scan frames"
          call fitsline(iufred,line)
          line = "BUNIT   = 'erg / s cm2 cm-1 sr'"
          call fitsline(iufred,line)
        elseif(itype==(-2)) then
          line = "EXTNAME = 'SCAN-NOISE'" &
            //"  / calculated noise in each scan frame"
          call fitsline(iufred,line)
          line = "BUNIT   = 'erg / s cm2 cm-1 sr'"
          call fitsline(iufred,line)
        elseif(itype==(-3)) then
          line = "EXTNAME = '2-D-NOISE'" &
            //"  / noise in a long-slit or camera frame"
          call fitsline(iufred,line)
          line = "BUNIT   = 'erg / s cm2 cm-1 sr'"
          call fitsline(iufred,line)
        else
          print '("*** Undefined extension type",i4)', itype
        endif
        if(naxis>1) then
          write(line,'("CTYPE2  = ''pixel along slit''")')
          call fitsline(iufred,line)
        endif
        write(line,'("END")')
        call fitsline(iufred,line)
        call fitslpad(iufred)

        return
      end


      subroutine tblxtend(iufitf,ns,nt,arr)
!	write the header and data for a fits table extension

        use ius
        use dims

        real arr(ms,nt)
        character(80) line
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf

        if(.not.dofits) return

        line = "XTENSION= 'TABLE   '"
        call fitsline(iufred,line)
        write(line,'("BITPIX  = ",i20)') 8
        call fitsline(iufred,line)
        write(line,'("NAXIS   = ",i20)') 2
        call fitsline(iufred,line)
!	table width
!        write(line,'("NAXIS1  = ",i20)') 12*nt
        write(line,'("NAXIS1  = ",i20)') 80
        call fitsline(iufred,line)
!	table length
        write(line,'("NAXIS2  = ",i20)') ns
        call fitsline(iufred,line)
        write(line,'("PCOUNT  = ",i20)') 0
        call fitsline(iufred,line)
        write(line,'("GCOUNT  = ",i20)') 1
        call fitsline(iufred,line)
!	entries in each line
        write(line,'("TFIELDS = ",i20)') nt
        call fitsline(iufred,line)
        line = "EXTNAME = 'EXTRACTED'"
        call fitsline(iufred,line)
        line = "TTYPE1  = 'WAVENUMBER'  / observed vacuum wavenumber"
        call fitsline(iufred,line)
        write(line,'("TBCOL1  = ",i20)') 1
        call fitsline(iufred,line)
        line = "TFORM1  = 'F12.4'"
        call fitsline(iufred,line)
        line = "TUNIT1  = 'cm^-1   '"
        call fitsline(iufred,line)
        line = "COMMENT = 'not corrected for motion of Earth or Sun'"
        call fitsline(iufred,line)
        line = "TTYPE2  = 'FLUX    '  / estimated flux through slit"
        call fitsline(iufred,line)
        write(line,'("TBCOL2  = ",i20)') 13
        call fitsline(iufred,line)
        line = "TFORM2  = 'E12.4'"
        call fitsline(iufred,line)
        line = "TUNIT2  = 'Jy      '"
        call fitsline(iufred,line)
        line = "TTYPE3  = 'NOISE   '  / estimated flux uncertainty"
        call fitsline(iufred,line)
        write(line,'("TBCOL3  = ",i20)') 25
        call fitsline(iufred,line)
        line = "TFORM3  = 'E12.4'"
        call fitsline(iufred,line)
        line = "TUNIT3  = 'Jy      '"
        call fitsline(iufred,line)
        line = "TTYPE4  = 'ATMO    '  / estimated atmo transmission"
        call fitsline(iufred,line)
        write(line,'("TBCOL4  = ",i20)') 37
        call fitsline(iufred,line)
        line = "TFORM4  = 'E12.4'"
        call fitsline(iufred,line)
        line = "COMMENT = 'from (blackbody-sky)/blackbody'"
        call fitsline(iufred,line)
        line = "TUNIT4  = 'none    '"
        call fitsline(iufred,line)
        line = "TTYPE5  = 'WAVELENGTH'  / observed vacuum wavelength"
        call fitsline(iufred,line)
        write(line,'("TBCOL5  = ",i20)') 49
        call fitsline(iufred,line)
        line = "TFORM5  = 'F12.6'"
        call fitsline(iufred,line)
        line = "TUNIT5  = 'micron  '"
        call fitsline(iufred,line)
        line = "END     "
        call fitsline(iufred,line)
        call fitslpad(iufred)

        do is = 1,ns
          write(line,'(f12.4,3es12.4,f12.6,20x)') (arr(is,ii),ii=1,5)
          call fitsline(iufred,line)
        enddo
        call fitslpad(iufred)

        return
      end


      subroutine fitsopen(iufitf,fitsfile,ierr)
!	open a fits data file or (in entries) write to an opened file

        use ius
        character(48) fitsfile
        character(80) line

        if(iufitf==iufraw) then
          call frawopen(iufitf,fitsfile,ierr)
        elseif(iufitf==iufred) then
          call fredopen(iufitf,fitsfile,ierr)
        else
          print '(" Invalid fits unit")'
          ierr = 9
        endif

        return

      entry fitsline(iufitf,line)

        if(iufitf==iufraw) then
          call frawline(iufitf,line)
        elseif(iufitf==iufred) then
          call fredline(iufitf,line)
        else
          print '(" Invalid fits unit")'
          ierr = 9
        endif

        return

      entry fitslpad(iufitf)

        if(iufitf==iufraw) then
          call frawlpad(iufitf)
        elseif(iufitf==iufred) then
          call fredlpad(iufitf)
        else
          print '(" Invalid fits unit")'
          ierr = 9
        endif

        return

      entry fitsclos(iufitf)

        if(iufitf==iufraw) then
          call frawclos(iufitf)
        elseif(iufitf==iufred) then
          call fredclos(iufitf)
        else
          print '(" Invalid fits unit")'
          ierr = 9
        endif

        return

      end


      subroutine frawopen(iufitf,fitsfile,ierr)
!	open a raw fits data file or (in entries) write to an opened file

        use ius
        use dims

        integer, parameter :: ITNONE = 0
        integer, parameter :: ITHEAD = 1
        integer, parameter :: ITDATA = 2

        character(48) fitsfile
        character(80) line,lines(mlines)
        logical leftopen
        real*4 rdata(mdata)
        real*4 rdat(1)
        integer*4 intdata(mdata),intdat(mdata)
        logical doflip,doflop,baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /byteflip/ doflip,doflop

        save nbytes,irec,itype
        save lines,intdata,rdata

        if(.not.(dofits.or.rdfits)) return
        if(index(fitsfile,' ')<=1) then
          ierr = 9
          return
        endif
        inquire(unit=iufitf,opened=leftopen)
        if(leftopen) then
          print '("*** Fits file",i3," was left open")', iufitf
          print '(" Last writing type",i2,2i6)', itype,irec,nbytes
          close(unit=iufitf)
        endif
        open(unit=iufitf,file=fitsfile,access='direct',recl=mbytes)
        nbytes = 0
        irec = 0
        itype = ITHEAD
        ierr = 0

        return

      entry frawline(iufitf,line)

        if(.not.dofits) return
        inquire(unit=iufitf,opened=leftopen)
        if(.not.leftopen) &
          print '(" Writing header to closed fits file",i3)', iufitf
        if(itype==ITNONE) then
          itype = ITHEAD
        elseif(itype>ITHEAD) then
          print '("*** Writing fits header line in data section")'
          itype = ITHEAD
        endif
        nbytes = nbytes+80
        iline = nbytes/80
        lines(iline) = line
        if(nbytes>=mbytes) then
          irec = irec+1
          write(unit=iufitf,rec=irec) (lines(i),i=1,mlines)
          nbytes = 0
        endif

        return

      entry frawlpad(iufitf)

        if(.not.dofits) return
        inquire(unit=iufitf,opened=leftopen)
        if(.not.leftopen) &
          print '(" Writing header to closed fits file",i3)', iufitf
        if(nbytes==0) then
          itype = ITNONE
          return
        elseif(itype/=ITHEAD) then
          print '("*** Attempt to pad fits header" &
            " when not writing to header")'
          return
        endif
        nlpad = (mbytes-nbytes)/80
        do ilpad = 1,nlpad
          nbytes = nbytes+80
          iline = nbytes/80
          lines(iline) = ' '
        enddo   
        irec = irec+1
        write(unit=iufitf,rec=irec) (lines(i),i=1,mlines)
        nbytes = 0
        itype = ITNONE

        return

      entry fitsidat(iufitf,intdat,ndat)

        if(.not.dofits) return
        inquire(unit=iufitf,opened=leftopen)
        if(.not.leftopen) &
          print '(" Writing ints to closed fits file",i3)', iufitf
        if(itype==ITNONE) then
          itype = ITDATA
        elseif(itype==ITHEAD) then
          print '("*** Writing fits data in header section")'
          itype = ITDATA
        endif
        idata = nbytes/4
        do idat = 1,ndat
          idata = idata+1
          nbytes = nbytes+4
          intdata(idata) = intdat(idat)
          if(nbytes>=mbytes) then
            irec = irec+1
            if(doflop) call flipi4(intdata,mdata)
            write(unit=iufitf,rec=irec) intdata
            idata = 0
            nbytes = 0
          endif
        enddo   
        
        return

      entry fitsipad(iufitf)

        if(.not.dofits) return
        inquire(unit=iufitf,opened=leftopen)
        if(.not.leftopen) &
          print '(" Writing ints to closed fits file",i3)', iufitf
        if(nbytes==0) then
          itype = ITNONE
          return
        elseif(itype/=ITDATA) then
          print '("*** Attempt to pad fits data section" &
            " when not writing data")'
          return
        endif
        idata1 = (nbytes/4)+1
        do idata = idata1,mdata
          intdata(idata) = 0
          nbytes = nbytes+4
        enddo   
        irec = irec+1
        if(doflop) call flipi4(intdata,mdata)
        write(unit=iufitf,rec=irec) intdata
        nbytes = 0
        itype = ITNONE

        return

      entry frawclos(iufitf)

        if(.not.dofits) return
        if((nbytes/=0).or.(itype/=ITNONE)) then
          print '("*** Closing unfinished fits file",2i4)',itype,nbytes
        endif
        close(unit=iufitf)

        return

      end



      subroutine fredopen(iufitf,fitsfile,ierr)
!	open a red fits data file or (in entries) write to an opened file

        use ius
        use dims

        integer, parameter :: ITNONE = 0
        integer, parameter :: ITHEAD = 1
        integer, parameter :: ITDATA = 2

        character(48) fitsfile
        character(80) line,lines(mlines)
        logical leftopen
        real*4 rdata(mdata)
        real*4 rdat(1)
        integer*4 intdata(mdata),intdat(mdata)
        logical doflip,doflop,baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /byteflip/ doflip,doflop

        save nbytes,irec,itype
        save lines,intdata,rdata

        if(.not.(dofits.or.rdfits)) return
        if(index(fitsfile,' ')<=1) then
          ierr = 9
          return
        endif
        inquire(unit=iufitf,opened=leftopen)
        if(leftopen) then
          print '("*** Fits file",i3," was left open")', iufitf
          print '(" Last writing type",i2,2i6)', itype,irec,nbytes
          close(unit=iufitf)
        endif
        open(unit=iufitf,file=fitsfile,access='direct',recl=mbytes)
        nbytes = 0
        irec = 0
        itype = ITHEAD
        ierr = 0

        return

      entry fredline(iufitf,line)

        if(.not.dofits) return
        inquire(unit=iufitf,opened=leftopen)
        if(.not.leftopen) &
          print '(" Writing header to closed fits file",i3)', iufitf
        if(itype==ITNONE) then
          itype = ITHEAD
        elseif(itype>ITHEAD) then
          print '("*** Writing fits header line in data section")'
          itype = ITHEAD
        endif
        nbytes = nbytes+80
        iline = nbytes/80
        lines(iline) = line
        if(nbytes>=mbytes) then
          irec = irec+1
          write(unit=iufitf,rec=irec) (lines(i),i=1,mlines)
          nbytes = 0
        endif

        return

      entry fredlpad(iufitf)

        if(.not.dofits) return
        inquire(unit=iufitf,opened=leftopen)
        if(.not.leftopen) &
          print '(" Writing header to closed fits file",i3)', iufitf
        if(nbytes==0) then
          itype = ITNONE
          return
        elseif(itype/=ITHEAD) then
          print '("*** Attempt to pad fits header" &
            " when not writing to header")'
          return
        endif
        nlpad = (mbytes-nbytes)/80
        do ilpad = 1,nlpad
          nbytes = nbytes+80
          iline = nbytes/80
          lines(iline) = ' '
        enddo   
        irec = irec+1
        write(unit=iufitf,rec=irec) (lines(i),i=1,mlines)
        nbytes = 0
        itype = ITNONE

        return

      entry fitsrdat(iufitf,rdat,ndat)

        if(.not.dofits) return
        inquire(unit=iufitf,opened=leftopen)
        if(.not.leftopen) &
          print '(" Writing reals to closed fits file",i3)', iufitf
        if(itype==ITNONE) then
          itype = ITDATA
        elseif(itype==ITHEAD) then
          print '("*** Writing fits data in header section")'
          itype = ITDATA
        endif
        idata = nbytes/4
        do idat = 1,ndat
          idata = idata+1
          nbytes = nbytes+4
          rdata(idata) = rdat(idat)
          if(nbytes>=mbytes) then
            irec = irec+1
            if(doflop) call flipr4(rdata,rdata,mdata)
            write(unit=iufitf,rec=irec) rdata
            idata = 0
            nbytes = 0
          endif
        enddo   
        
        return

      entry fitsrpad(iufitf)

        if(.not.dofits) return
        inquire(unit=iufitf,opened=leftopen)
        if(.not.leftopen) &
          print '(" Writing reals to closed fits file",i3)', iufitf
        if(nbytes==0) then
          itype = ITNONE
          return
        elseif(itype/=ITDATA) then
          print '("*** Attempt to pad fits data section" &
            " when not writing data")'
          return
        endif
        idata1 = (nbytes/4)+1
        do idata = idata1,mdata
          rdata(idata) = 0.
          nbytes = nbytes+4
        enddo   
        irec = irec+1
        if(doflop) call flipr4(rdata,rdata,mdata)
        write(unit=iufitf,rec=irec) rdata
        nbytes = 0
        itype = ITNONE

        return

      entry fredclos(iufitf)

        if(.not.dofits) return
        if((nbytes/=0).or.(itype/=ITNONE)) then
          print '("*** Closing unfinished fits file",2i4)',itype,nbytes
        endif
        close(unit=iufitf)

        return

      end



      subroutine fitsarr(arr,fitsfile,ierr)
!       store reduced camera, long-slit, or untorted data array to fitsfile

        use ius
        use dims

        real arr(mx,my)
        character(48) fitsfile
        character(60) bunits,comment

        common /nn/ nx,ny,nc,norder,ns,nt

        ierr = 0
        bunits = 'erg/s cm2 sr cm-1'
        comment = 'cgs intensity units'
        call fithchar('BUNIT   ',bunits,comment,iufith)
        if(index(fitsfile,' ')<=1) return
        call redstart(nx,ny-4,1,fitsfile,iufith,ierr)
        if(ierr>0) return
 
        do iy = 5,ny
          call fitsrdat(iufred,arr(1,iy),nx)
        enddo
        call fitsrpad(iufred)

        do iy = 1,4
          call redxtend(iufred,nx,1,0,iy)
          call fitsrdat(iufred,arr(1,iy),nx)
          call fitsrpad(iufred)
        enddo

        call fitsclos(iufred)
        return

      end



      subroutine fitspec(spec,fitsfile,ierr)
!       store xd spec array as fits
!	with first four columns stored in extensions

        use ius
        use dims
        use modes

        real spec(ms,mt)
        character(48) fitsfile
        character(60) comment,bunits

        common /nn/ nx,ny,nc,norder,ns,nt

        ierr = 0
        bunits = 'erg/s cm2 sr cm-1'
        comment = 'cgs intensity units'
        call fithchar('BUNIT   ',bunits,comment,13)
        if(index(fitsfile,' ')<=1) return
        call redstart(ns,nt-4,1,fitsfile,iufith,ierr)
        if(ierr>0) return

        do it = 5,nt
          call fitsrdat(iufred,spec(1,it),ns)
        enddo
        call fitsrpad(iufred)

        do it = 1,4
          call redxtend(iufred,ns,1,0,it)
          call fitsrdat(iufred,spec(1,it),ns)
          call fitsrpad(iufred)
        enddo   

        call fitsclos(iufred)
        return

      end



      subroutine fitspeci(speci,fitsfile,ierr)
!       store ny=4 spec array as fits (for drain)

        use ius
        use dims
        use modes

        real speci(ms,4)
        character(48) fitsfile

        common /nn/ nx,ny,nc,norder,ns,nt

        ierr = 0
        if(index(fitsfile,' ')<=1) return
        call redstart(1,1,1,fitsfile,iufith,ierr)
        if(ierr>0) return

        do it = 1,4
          call redxtend(iufred,ns,1,0,it)
          call fitsrdat(iufred,speci(1,it),ns)
          call fitsrpad(iufred)
        enddo

        call fitsclos(iufred)
        return

      end



      subroutine fitscan(scan,ifsum,fitsfile,ierr)
!       store a scan array as a 3-D fits array in ny x nz x nx order
!	with iy=1-4 stored in extensions

        use ius
        use dims

        real scan(mu,mv,mw)
        real temp(mv)
        character(48) fitsfile
        character(60) comment,bunits
        logical ifsum

        common /nn/ nx,ny,nc,norder,ns,nt

        ierr = 0
        if(ifsum) then
          lv = nvz-4
          lw = nwz
        else
          lv = nv-4
          lw = nw
        endif
        bunits = 'erg/s cm2 sr cm-1'
        comment = 'cgs intensity units'
        call fithchar('BUNIT   ',bunits,comment,13)
        if(index(fitsfile,' ')<=1) then
          ierr = 1
          return
        endif
        if(ifsum) then
          call redstart(nu,lv,lw,fitsfile,iufish,ierr)
        else
          call redstart(nu,lv,lw,fitsfile,iufith,ierr)
        endif
        if(ierr>0) return

        do iu = 1,nu
          do iw = 1,lw
            do iv = 1,lv
              temp(iv) = scan(iu,iv+4,iw)
            enddo
            call fitsrdat(iufred,temp,lv)
          enddo
        enddo
        call fitsrpad(iufred)

        do iv = 1,4
          call redxtend(iufred,nu,1,0,iv)
          call fitsrdat(iufred,scan(1,iv,1),nu)
          call fitsrpad(iufred)
        enddo

        call fitsclos(iufred)
        return

      end


      subroutine fitstotx(fitsfile,txfile,ierr)

        use ius

        character(48) fitsfile
        character(32) txfile

        call fitsopen(iufraw,fitsfile,ierr)
        if(ierr/=0) then
          print '(" Error opening ",a48)', fitsfile
          return
        endif
        call fitsrdhd(iufraw,irec,ierr)
        if(ierr/=0) then
          print '(" Error reading header ",a48)', fitsfile
          close(unit=iufraw)
          return
        endif
        call fitsread(iufraw,irec,txfile,ierr)
        if(ierr/=0) then
          print '(" Error reading data ",a48)', fitsfile
          close(unit=iufraw)
          return
        endif
        close(unit=iufraw)

      end
