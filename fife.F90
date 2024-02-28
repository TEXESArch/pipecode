#include "mods.h"

      program fife

!	program to read in raw texes data,
!	process, and write back to disk
!	also write raw and reduced data in fits format

!     subroutine readscript(kbd,pipedir,fitsdir,ierr)
!	read commands either from pipescript or from the keyboard
!     subroutine settorts(date)
!	initialize some tort parameters
!     subroutine waitasec(pause,image,arr,nx,lx,ny)
!	if(pause) let the user find a value of a pixel
!     subroutine cook(fname,seqno,pipedir,ierr)
!	the main data reduction routine
!     subroutine makenames(fname,seqno,pipedir,fitsdir)
!	make raw and reduced file names
!     subroutine setgood(date,ierr)
!	initialize good array, flagging known bad pixels
!     subroutine setone(arr,nz,ierr)
!	set an array to 1.
!     subroutine readhdr(ifobj,ierr)
!	read the raw header and copy to red header and fits.hd
!     subroutine vels(day,ra,dec)
!	calculate Earth's and Sun's motions
!     subroutine precess(ra,dec,y1,y2)
!	precess coordinates for vels
!     subroutine readcard(cardfile,ierr)
!	read a flat file
!     subroutine readraw(rawfile,abeam,bbeam,np,ierr)
!	read a raw data file
!     subroutine makeflat(flat,ierr)
!	make a flat field from a flat file; also run testtort
!     subroutine oldspike(arr,var,varfac,nz,ierr)
!	old routine to look for clouds and spikes
!     subroutine despike(arr,var,varfac,nz,ierr)
!	must be better than the old version
!     subroutine debounce(arr,brr,flat,nz,ierr)
!	try to remove 1st and 2nd derivative bounce by shifting and smoothing
!     subroutine diffarr(a,b,d,avar,bvar,std,nz,ierr)
!	subtract A-B beams and take noise from quieter
!     subroutine calibrate(arr,flat,nz,ierr)
!	multiply an array by flat to calibrate array
!     subroutine clean(arr,std,nz,ierr)
!	interpolate over bad and noisy pixels
!     subroutine tort(a,nz,ierr)
!	remove optical distortions
!     subroutine qtort(a,nz,ierr)
!	a quick tort for debounce
!     subroutine testtort(c,ierr)
!	try out tort on black-shiny to find parameters
!     subroutine submean(arr,flat,nz,ierr)
!	make slit mean zero at each spectral point
!     subroutine subcorr(arr,brr,flat,nz,ierr)
!	for map mode, remove correlated sky noise
!     subroutine cirrus(arr,brr,flat,nz,ierr)
!	for nod-off mode, subtract sky mean and slope along orders
!     subroutine setillum(array,ierr)
!	find which pixels are on array and in orders
!     subroutine newillum(array,ierr)
!	find which pixels are on array and in orders (new version)
!     subroutine maketempl(arr,flat,std,templ,nz,
!	make template for shift, wtadd, and extract
!     subroutine shift(arr,flat,std,nz,ierr)
!	shift diffs spatially to match sum or sumspec
!     subroutine wtadd(arr,sumarr,flat,std,nz,ierr)
!	weight diffs by correlation with straight sum, and average
!     subroutine adddiffs(arr,sumarr,std,nz,ierr)
!	average diffs with equal weighting
!     subroutine extract(arr,flat,std,card,spec,ierr)
!	extract 1-d spectrum from 2-d data; also calculate wnos
!     subroutine zerosum(ierr)
!	initialize sumspec and scansum to zero
!     subroutine sumarr(arr,ierr)
!	add long-slit spectrum to arrsum
!     subroutine sumspec(spec,ierr)
!	add xd spectrum to specsum
!     subroutine sumscan(scan,scansum,std,nz,im,sumtime,ierr)
!	add scan to scansum, with optional shift
!     subroutine subscansky(scan,std,nz,ierr)
!	subtract interpolated sky from scan
!     subroutine subscancorr(scan,std,nz,ierr)
!	subtract correlated sky fluctuations from scan if(skyint(1)<0)
!     subroutine scanmap(scan,std,nz)
!	make postage stamp scan images
!     subroutine writeparms(iunit,junit,ierr)
!	write parameters to reduced header
!     subroutine storearr(arr,std,outfile,fitsfile,ierr)
!	store untorted or long-slit array
!     subroutine storescan(scan,flat,std,nz,outfile,fitsfile,
!	sumtime,nssum,ierr)
!	store scan data
!     subroutine storespec(spec,outfile,fitsfile,ierr)
!	store xd spectrum
!     subroutine storesum(outfile,fitsfile,ierr)
!	store summed xd or ls spectrum
!     subroutine showsum(ierr)
!	plot summed spectrum
!     subroutine smooth(nx,ny,arr,sarr,sx,sy)
!	seems to be unused
!     subroutine checkarr(arr,name,ierr)
!	check for NaN or Inf in arr
!     subroutine matinv(arr,ni,det)
!	invert a symmetric array

        use ius
        use dims
        use modes
        use paths

        logical good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        character(60) object,feature,obsmode,instmode
        character(60) objtype,piname,pid,note,weather,warning,version
        real nodpa,lores,kmirror,krot,lskrot
        character(32) rawfile,cardfile,lastcard,redfile,sumfile,posfile
        character(48) rawfits,cardfits,redfits,sumfits
        character(4) ired
        character(8) pipedir,fitsdir
        character(32) pipelog
        character(40) wnoline
        character(8) yn
        logical kbd,doplot
        integer pgopen
        logical doarch,scriptopen
        logical doflip,doflop

        common /byteflip/ doflip,doflop
        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
               piname,pid,note,weather,warning,version
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
        common /tort0/ hrfl0,xdfl0,xdmrdgr,xdlrdgr,xdmr0,xdlr0,fred0, &
                xdkrot,lskrot
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /filenames/rawfile,cardfile,lastcard,redfile,sumfile,posfile, &
                rawfits,cardfits,redfits,sumfits,ired
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime
        common /iwin/ iwin1,iwin2

        version = '23.0904'
        print '(" Version ",a8," of fife(c) has been certified bug-free")', &
          version
        print '(" And conforming to ASTM A53 pipe standards")'

!	flip on input (with intel inside(c))
        doflip = .true.
!	flip on output (your choice)
        doflop = doflip

        nx = mx
        ny = my
        nc = 3

!	default header values
        pid = ' '
        piname = ' '
        weather = ' '
        warning = ' '
        waveno0 = 1000.
        temp = 273.16
        airmass = 1.
        frtime = 1.
        obstime = 0.
        tottime = 0.
        beamtime = 0.
        addtime = 0.
        echelle = 60.
        lores = 20.
        gain = 1.0
        pixelwd = 0.0030
        slit = 180.
        efl = 12.*300./2.
        nbitpix = 32

!	default reduction parameters
        modecard = MBLKSKY
        modeblk = MUNKNOWN
        modetort = MBLK
        modeext = MUNWT
        intext(1) = 0
        intext(2) = 0
        intext(3) = 0
        intext(4) = 0
        lastcard = ' '
        doplot = .true.
        verbose = .false.
        ask = .false.
        pause = .false.
        sincwt = .true.
        doshift = .false.
        intshift(1) = 0
        intshift(2) = 0
        intshift(3) = 0
        intshift(4) = 0
        intshift(5) = 0
        intshift(6) = 0
        doaddwt = .false.
        dosubsky = .false.
        intsky(1) = 0
        intsky(2) = 0
        intsky(3) = 0
        intsky(4) = 0
        intill(1) = 0
        intill(2) = 0
        intill(3) = 0
        intill(4) = 0
        dosum = .false.
        doscansum = .false.
        abba = .false.
        crossdisp = .true.
        ffttort = .false.
        kfix = .false.
        baddata = .false.
        ierr = 0
        call zerosum(ierr)
        thrfac = 0.5
        spikefac = 30.
        stdfac = 20.
        satval = 0.94*2.**16
        xnlin = 0.
        ynlin = 0.
        znlin = 0.
        dtcal = 0.
        dtsky = 0.
        darkval = 0.
        cloud = 0.
        bounce = 0.
        slitpa = 999.
        objtype = 'unkn'
        piname = ' '
        pid = ' '
        do iz = 1,mp
          wt(iz) = 1.
        enddo
        ired = ' '
        fitsdir = ' '
        posfile = ' '
        domon = .false.
        dofits = .false.
        rdfits = .false.

!	default distortion parameters (some overridden when date set)
        date = 0.0
        slitrot = 0.0
        krot = 0.0
        lskrot = 0.0
        xdkrot = 0.0
        detrot = 0.080
        hrr = 9.8
        hrg = 0.010
        hrdgr = 0.7590
        hrdgr = 0.3*2.54*0.996*cos(hrg)
        xdg = 0.0133
!	xddgr and fls set in readheader from instmode and waveno0
!	did I have these set small to compensate for an error in xdfl?
!	xdmrdgr = 0.003143
        xdmrdgr = 0.003151
!	xdlrdgr = 0.001318
        xdlrdgr = 0.001328
!	I don''t think the displacement of the xd parab changes efl.
!	hrfl0 was 40 inches.  I think it should be shrunken by .996
!	note that it may get overwritten in settorts
!	and there is some evidence that the fls differ from these values
!	xdfl0 should be a shrunken 30 inches
        hrfl0 = 101.2
        xdfl0 = 75.9
        xddgr = xdmrdgr
!	hrfl and xdfl will be recalculated based on wavelength
        hrfl = hrfl0/2.1
        xdmr0 = 0.2  !  Why?
        xdlr0 = 0.2
        xdfl = xdfl0/2.1
        fred0 = 2.12
        spacing = 32.
        norder = 1
        ns = 256
        nt = mt
        wno0 = 1.0
        radvel = 0.0
        call veltype(0)
        brl = 1.5
        x0brl = 50.
        y0brl = 0.

!	open script and log files

        print '(" Enable plotting? (y) ",$)'
        read '(a8)', yn
        doplot = (index(yn,'n')==0)
        print '(" Use keyboard command mode? (n) ",$)'
        read '(a8)', yn
        kbd = (index(yn,'y')>0)
        if(.not.kbd) then
          inquire(file='pipescript_archive',exist=doarch)
          if(doarch) then
            open(unit=iups,file='pipescript_archive',iostat=ierr)
            print '(" Taking commands from pipescript_archive")'
          else
            inquire(file='archscript',exist=doarch)
            if(doarch) then
              open(unit=iups,file='archscript',iostat=ierr)
              print '(" Taking commands from archscript")'
            else
              open(unit=iups,file='pipescript',iostat=ierr)
              if(ierr==0) then
                print '(" Taking commands from pipescript")'
              else
                print '(" Error opening pipescript")'
                kbd = .true.
              endif
            endif
          endif
        endif
        pipedir = 'pipe/'
        print '(" Store to directory pipe?  (y) ",$)'
        read '(a8)', yn
        if(index(yn,'n')>0) then
          do
            print '(" Enter directory name ",$)'
            read '(a8)', pipedir
            last = index(pipedir,' ')-1
            if((last>0).and.(last<8)) exit
            print '(" Directory name must be 1-7 characters long")'
          enddo
          last = min(last,7)
          pipedir = pipedir(1:last)//'/'
        endif
        last = index(pipedir,' ')-1
        pipelog = pipedir(1:last)//'pipelog'
        print '(" Overwrite pipelog? (y) ",$)'
        read '(a8)', yn
        if(index(yn,'n')>0) then
          open(unit=iupl,file=pipelog,access='append',iostat=ierr)
        else
          open(unit=iupl,file=pipelog,iostat=ierr)
        endif
        if(ierr/=0) then
          print '(" pipe directory not found.", &
            "  Putting pipelog in current directory.")'
          if(index(yn,'n')>0) then
            open(unit=iupl,file='pipelog',access='append',iostat=ierr)
          else
            open(unit=iupl,file='pipelog',iostat=ierr)
          endif
        endif
        open(unit=iuwno,file='pipewnos')

!	start pgplot

        if(doplot) then
!	pgdev should be set in mods.h
          iwin1 = pgopen(pgdev)
          call pgask(.false.)
          if(iwin1<=0) then
            print '(" Error opening xwindow 1 for pgplot")'
            print '(" Try changing pgopen to /xserve or /xwindow")'
            print '(" And while you''re at it check atmodir")'
            go to 990
          endif
          call pgpap(7.5,1.0)
          call pgenv(0.,1.,0.,1.,1,-2)
          iwin2 = iwin1
          print '(" Arrange your windows and hit RETURN ",$)'
          read '(a8)', yn
        else
          iwin1 = -1
          iwin2 = -1
        endif

!	read script file or keyboard for further instructions

        do
          call readscript(kbd,pipedir,fitsdir,ierr)
          if(ierr<0) then
            exit
          elseif(ierr>0) then
            if(kbd) then
              print '(" Error in command line")'
            else
              print '(" Error in readscript.  ",i2)', ierr
              kbd = .true.
            endif
          endif
        enddo
        close(unit=iuwno)
        open(unit=iuwno,file='pipewnos')
        print '(" Fitted wno0,hrr values")'
        do
          read(unit=iuwno,fmt='(a40)',iostat=ierr) wnoline
          if(ierr>0) exit
          print '(a40)', wnoline
        enddo

        call pgend
  990   close(unit=iups)
        close(unit=iupl)
        close(unit=iuwno)

        stop
      end


      subroutine quit
        use ius

        call pgend
        close(unit=iups)
        close(unit=iupl)
        close(unit=iuwno)

        stop
      end


      subroutine readscript(kbd,pipedir,fitsdir,ierr)
!	read a script or kbd line and do what it says

        use ius
        use dims
        use modes

        logical kbd
        character(8) pipedir,fitsdir
        character(80) line,pname,pval
        character(24) fname,fnold,gotoline
        character(4) seqno,ired
        character(32) rawfile,cardfile,lastcard,redfile,sumfile,posfile
        character(48) rawfits,cardfits,redfits,sumfits
        character(60) object,feature,obsmode,instmode
        character(60) objtype,piname,pid,note,weather,warning,version
        character(24) scriptfile
        character(8) yn
        real nodpa,lores,kmirror,krot,lskrot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical good
        logical scriptopen,sumopen,doarch
        logical logstr
        logical doflip,doflop
        logical :: checkwno = .true.	! check wno0 for each new object
        logical :: checkall = .false.	! check wno0 for every file

        common /byteflip/ doflip,doflop
        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
                piname,pid,note,weather,warning,version
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
        common /tort0/ hrfl0,xdfl0,xdmrdgr,xdlrdgr,xdmr0,xdlr0,fred0, &
                xdkrot,lskrot
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /filenames/rawfile,cardfile,lastcard,redfile,sumfile,posfile, &
                rawfits,cardfits,redfits,sumfits,ired
        common /iwin/ iwin1,iwin2
        save checkwno,checkall

        ierr = 0
        if(date==0.) crossdisp = .true.
        if(kbd) then
          print '(" fife > ",$)'
          read '(a79)', line
        else
          read(unit=iups,fmt='(a79)',iostat=ierr) line
          if(ierr==0) then
            print '(a80)', line
          elseif(ierr<0) then
            print '(" Read to end of script")'
            kbd = .true.
            ierr = 9
            return
          elseif(ierr>0) then
            print '(" Error reading script")'
            ierr = 9
            return
          endif
        endif
        write(iupl,'(a80)') line
        last = 78
        do jrst = 1,last
          irst = jrst
          if(index(line(irst:last),'#')==1) return
          if(index(line(irst:last),' ')==1) cycle
          if(index(line(irst:last),'	')==1) cycle ! that must be tab
          exit
        enddo
        if(irst>=last) return
        line = line(irst:last)
        if(index(line,'#')>0) then
          last = index(line,'#')-1
          if(last<1) return
          line = line(irst:last)
        endif
        if(index(line,'goto')>0) then
          irst = index(line,'goto')+4
          do jrst = irst,last
            if(index(line(jrst:last),' ')==1) cycle
            irst = jrst
            exit
          enddo
          do jast = irst,last
            if(index(line(jast:last),' ')==1) then
              last = jast-1
              exit
            endif
          enddo
          gotoline = line(irst:last)
          last = last-irst+1
          if(kbd) then
            inquire(unit=iups,opened=scriptopen)
            if(.not.scriptopen) then
              print '(" Enter scriptfile name: "$)'
              read '(a24)', scriptfile
              open(unit=iups,file=scriptfile)
            endif
            kbd = .false.
          endif
          do
            read(unit=iups,fmt='(a79)',iostat=ierr) line
            if(ierr/=0) then
              print '(" Error or end of script.  Can''t read to ",a24)', &
                gotoline
              kbd = .true.
              ierr = 9
              return
            endif
            if(index(line,gotoline(1:last))>0) exit
          enddo
          return
        endif
        ierr = 0
        if(index(line,'end')==1) then
          ierr = -1
          if(kbd) then
            print '(" Hasta luego")'
          else
            print '(" End of script")'
          endif
          return
        elseif((index(line,'bye')==1).or.(index(line,'quit')==1) &
          .or.(index(line,'exit')==1)) then
          ierr = -1
          print '(" Adios")'
          return
        elseif(index(line,'buy')==1) then
          ierr = -1
          print '(" Sell")'
          return
        elseif(index(line,'kbd')==1) then
          print '(" Entering keyboard command mode")'
          kbd = .true.
          return
        elseif(index(line,'cont')==1) then
          print '(" Taking commands from script")'
          kbd = .false.
          return
        elseif(index(line,'script')==1) then
          ival = index(line,'=')+1
          if(ival>1) then
            do jval = ival,80
              if(jval==80) then
                scriptfile = 'pipescript'
                exit
              elseif(line(jval:jval)==' ') then
                cycle
              else
                kval = index(line(jval:79),' ')
                if(kval==0) kval = 80
                scriptfile = line(jval:jval+kval-1)
                exit
              endif
            enddo
          else
            inquire(file='pipescript_archive',exist=doarch)
            if(doarch) then
              scriptfile = 'pipescript_archive'
            else
              inquire(file='archscript',exist=doarch)
              if(doarch) then
                scriptfile = 'archscript'
              else
                scriptfile = 'pipescript'
              endif
            endif
          endif
          print '(" Taking commands from ",a24)', scriptfile
          kbd = .false.
          inquire(unit=iups,opened=scriptopen)
          if(scriptopen) close(unit=iups)
          open(unit=iups,file=scriptfile)
          return
        elseif(index(line,'storesum')==1) then
          if((modeobs<MSCAN).and.(.not.testrun)) then
            call storesum(sumfile,sumfits,ierr)
          endif
          close(unit=iusumh)
          if(dofits) close(unit=iufish)
          sumtime = 0.
          dosum = .false.
        elseif(index(line,'showsum')==1) then
          if(modeobs<MSCAN) call showsum(ierr)
        elseif(index(line,'zerosum')==1) then
          call zerosum(ierr)
        elseif(index(line,'setwno')==1) then
          wno0 = 0.
          norder = 0
          kfix = .false.
        elseif(index(line,'longslit')==1) then
          if(date==0.) then
            print '("***That won''t work!  Read date or a file first.")'
            kbd = .true.
            return
          endif
          if(crossdisp) then
            xdkrot = krot
            if(lskrot/=0.) then
!	will be changed in lo-res mode
              krot = lskrot
            else
              print '("***Why hasn''t lskrot been set?")'
!              lskrot = .012
!	in case we really want lskrot = 0
              krot = lskrot
            endif
            if(intext(1)/=0) print &
              '("***WARNING extint was set in crossdisp mode")'
            crossdisp = .false.
            hrr = 9.8
          else
            print '("***Weren''t you already in longslit mode?")'
          endif
          detrot = .083
        elseif(index(line,'crossdisp')==1) then
          if(.not.crossdisp) then
!	would be wrong in lo-res mode
!            lskrot = krot
            if(xdkrot/=0.) then
              krot = xdkrot
            else
              print '("***Why hasn''t xdkrot been set?")'
              xdkrot = -.020
              krot = -.020
            endif
            if(intext(1)/=0) print &
              '("***WARNING extint was set in longslit mode")'
            crossdisp = .true.
            hrr = 9.8
          else
            print '("***Weren''t you already in crossdisp mode?")'
          endif
          detrot = .080
        elseif(index(line,'baddata')>0) then
          baddata = .true.
          print '(" flagging next data file as bad")'
        elseif(index(line,'=')>0) then
          iname = index(line,'=')-1
          pname = line(1:iname)
          ival = iname+2
    8     pval = line(ival:79)
          if(index(pval,' ')==1) then
            ival = ival+1
            if(ival<79) goto 8
            print '(" No value given")'
          endif
          ierr = 0
          if(index(pname,'PID')>0) then
            ival = index(pval,' ')-1
            if((ival<1).or.(ival>16)) then
              print '(" PID must be 1-16 characters long")'
              print '(" Enter PID")'
              read '(a16)', pval
            endif
            ival = index(pval,' ')-1
            ival = min(ival,16)
            if(ival>0) then
              pid = pval(1:ival)
            else
              pid = ' '
            endif
          elseif(index(pname,'PI')>0) then
            ival = index(pval,' ')-1
            if((ival<1).or.(ival>16)) then
              print '(" PI name must be 1-16 characters long")'
              print '(" Enter PI name")'
              read '(a16)', pval
            endif
            ival = index(pval,' ')-1
            ival = min(ival,16)
            if(ival>0) then
              piname = pval(1:ival)
            else
              piname = ' '
            endif
          elseif(index(pname,'note')>0) then
            note = pval
          elseif(index(pname,'weather')>0) then
            weather = pval
          elseif(index(pname,'warning')>0) then
            warning = pval
          elseif(index(pname,'flip')>0) then
            doflip = logstr(pval,ierr)
          elseif(index(pname,'flop')>0) then
            doflop = logstr(pval,ierr)
          elseif(index(pname,'newdir')>0) then
            call chdir(pval)
          elseif(index(pname,'pipedir')>0) then
            ival = index(pval,' ')-1
            if((ival<1).or.(ival>7)) then
              print '(" Directory name must be 1-7 characters long")'
              print '(" Enter directory name")'
              read '(a7)', pval
            endif
            ival = index(pval,' ')-1
            ival = min(ival,7)
            if(ival>0) then
              pipedir = pval(1:ival)//'/'
            else
              print '(" What?")'
            endif
          elseif(index(pname,'mon')>0) then
            domon = logstr(pval,ierr)
          elseif(index(pname,'fitsdir')>0) then
            ival = index(pval,' ')-1
            if((ival<1).or.(ival>7)) then
              print '(" Directory name must be 1-7 characters long")'
              print '(" or ''none''.  Enter directory name")'
              read '(a7)', pval
            endif
            ival = index(pval,' ')-1
            ival = min(ival,7)
            if((ival>0).and.(index(pval,'none')==0)) then
              fitsdir = pval(1:ival)//'/'
              dofits = .true.
            else
              print '(" Turning off fits saving")'
              dofits = .false.
            endif
          elseif(index(pname,'dofits')>0) then
            if((index(pval,'f')>0).or.(index(pval,'none')>0)) then
              dofits = .false.
            elseif(index(pval,'t')>0) then
              if(index(fitsdir,' ')<=1) then
                print '("*** fitsdir not initialized")'
                kbd = .true.
              else
                dofits = .true.
              endif
            else
              print '("*** What?")'
            endif
          elseif((index(pname,'rdfits')>0).or.(index(pname,'readf')>0)) then
            if(index(pval,'f')>0) then
              rdfits = .false.
            elseif(index(pval,'t')>0) then
              rdfits = .true.
            else
              print '("*** What?")'
            endif
          elseif(index(pname,'date')>0) then
            read(unit=pval,fmt=*,iostat=ierr) date
            if(ierr==0) then
              if(date>2000.) date = date-2000.
              iyr = int(date)
              imo = int(100.*(date-iyr))
              idy = nint(100.*(100.*(date-iyr)-imo))
              call settorts(date)
            endif
          elseif(index(pname,'posfile')>0) then
            ival = index(pval,' ')-1
            ival = min(ival,31)
            if((ival>0).and.(index(pval,'none')==0)) &
              posfile = pval(1:ival)
          elseif(index(pname,'slitpa')>0) then
            read(unit=pval,fmt=*,iostat=ierr) slitpa
            if(slitpa<0.) slitpa = slitpa+360.
            if(slitpa<720.) slitpa = -slitpa
          elseif(index(pname,'checkall')>0) then
            checkall = logstr(pval,ierr)
          elseif(index(pname,'checkwno')>0) then
            checkwno = logstr(pval,ierr)
          elseif((index(pname,'cardmode')>0) &
              .or.(index(pname,'flatmode')>0)) then
            if(index(pval,'blksk')>0) then
              modecard = MBLKSKY
            elseif(index(pval,'blkob')>0) then
              modecard = MBLKOBJ
            elseif(index(pval,'blksh')>0) then
              modecard = MBLKSHINY
            elseif(index(pval,'bl')>0) then
              modecard = MBLK
            elseif(index(pval,'sk')>0) then
              modecard = MSKY
            elseif(index(pval,'sh')>0) then
              modecard = MSHINY
            elseif(index(pval,'none')>0) then
              modecard = MNONE
              if((modetort>1).and.(modetort<8)) then
                print '(" Setting tortmode = obj")'
                modetort = MOBJ
              endif
            elseif(index(pval,'obj')>0) then
              modecard = MOBJ
            elseif((index(pval,'last')>0) &
              .or.(index(pval,'old')>0)) then
!	note: can''t just use old flat since it is torted
              if(.not.(testrun.or.rdfits)) modecard = -modecard
            elseif(index(pval,'bsbs')>0) then
              modecard = MBSBS
            else
              print '("*** Unknown flatmode: ",a)', pval
              modecard = MUNKNOWN
            endif
            if((abs(modecard)>4).and.(abs(modecard)<=8)) then
              modetort = MBLK
            elseif(modetort>1) then
              modetort = abs(modecard)
            endif
          elseif((index(pname,'blk')>0) &
            .or.(index(pname,'black')>0)) then
            read(unit=pval,fmt='(i8)',iostat=ierr) modeblk
          elseif(index(pname,'tortmode')>0) then
            if((index(pval,'card')>0) &
              .or.(index(pval,'flat')>0)) then
              modetort = modecard
            elseif(index(pval,'none')>0) then
              modetort = MNONE
            elseif(index(pval,'bl')>0) then
              modetort = MBLK
            elseif(index(pval,'sk')>0) then
              modetort = MSKY
            elseif(index(pval,'sh')>0) then
              modetort = MSHINY
            elseif(index(pval,'obj')>0) then
              modetort = MOBJ
            elseif(index(pval,'cell')>0) then
              modetort = MCELL
            elseif(index(pval,'old')>0) then
              modetort = MOLD
            else
              print '("*** Unknown tortmode: ",a)', pval
              modetort = MUNKNOWN
            endif
          elseif(index(pname,'thr')>0) then
            read(unit=pval,fmt=*,iostat=ierr) thrfac
          elseif(index(pname,'spike')>0) then
            read(unit=pval,fmt=*,iostat=ierr) spikefac
          elseif((index(pname,'stdfac')>0) &
              .or.(index(pname,'noise')>0)) then
            read(unit=pval,fmt=*,iostat=ierr) stdfac
          elseif(index(pname,'sat')>0) then
            read(unit=pval,fmt=*,iostat=ierr) satval
            satval = satval*2.**16
          elseif(index(pname,'dark')>0) then
            read(unit=pval,fmt=*,iostat=ierr) darkval
          elseif(index(pname,'xnlin')>0) then
            read(unit=pval,fmt=*,iostat=ierr) xnlin
          elseif(index(pname,'ynlin')>0) then
            read(unit=pval,fmt=*,iostat=ierr) ynlin
          elseif(index(pname,'znlin')>0) then
            read(unit=pval,fmt=*,iostat=ierr) znlin
          elseif(index(pname,'nlin')>0) then
            read(unit=pval,fmt=*,iostat=ierr) xnlin
          elseif(index(pname,'dtcal')>0) then
            read(unit=pval,fmt=*,iostat=ierr) dtcal
          elseif(index(pname,'dtblk')>0) then
            read(unit=pval,fmt=*,iostat=ierr) dtcal
          elseif(index(pname,'dtsky')>0) then
            read(unit=pval,fmt=*,iostat=ierr) dtsky
          elseif(index(pname,'cloud')>0) then
            read(unit=pval,fmt=*,iostat=ierr) cloud
          elseif(index(pname,'bounce')>0) then
            read(unit=pval,fmt=*,iostat=ierr) bounce
          elseif(index(pname,'slitrot')>0) then
            read(unit=pval,fmt=*,iostat=ierr) slitrot
          elseif(index(pname,'lskrot')>0) then
            read(unit=pval,fmt=*,iostat=ierr) lskrot
            print '(" Setting mrkrot =",f6.3," lrkrot =",f6.3)', &
              lskrot,lskrot-0.004
!	will be changed in lo-res mode
            if(.not.crossdisp) krot = lskrot
          elseif(index(pname,'xdkrot')>0) then
            read(unit=pval,fmt=*,iostat=ierr) xdkrot
            if(crossdisp) krot = xdkrot
          elseif(index(pname,'krot')>0) then
            read(unit=pval,fmt=*,iostat=ierr) krot
            if(crossdisp) then
              xdkrot = krot
            else
              lskrot = krot
              print '(" Setting mrkrot =",f6.3," lrkrot =",f6.3)', &
                lskrot,lskrot-0.004
            endif
          elseif(index(pname,'detrot')>0) then
            read(unit=pval,fmt=*,iostat=ierr) detrot
          elseif(index(pname,'hrfl')>0) then
            read(unit=pval,fmt=*,iostat=ierr) hrfl0
          elseif(index(pname,'fred')>0) then
            read(unit=pval,fmt=*,iostat=ierr) fred0
          elseif(index(pname,'hrr')>0) then
            read(unit=pval,fmt=*,iostat=ierr) hrr
          elseif(index(pname,'hrg')>0) then
            read(unit=pval,fmt=*,iostat=ierr) hrg
          elseif(index(pname,'hrdgr')>0) then
            read(unit=pval,fmt=*,iostat=ierr) hrdgr
          elseif(index(pname,'wno')>0) then
            read(unit=pval,fmt=*,iostat=ierr) wno0
            if(modecard==MNONE) then
              wno0 = abs(wno0)
            else
              wno0 = -wno0
            endif
            kfix = .false.
            norder = 0
            hrr0 = hrr
            wnoa = abs(wno0)
!	crossdisp may change when a file is read, but this shouldn't hurt
            if(crossdisp.and.(wnoa*hrr/=0.)) then
              hrorder = nint(2.0*hrdgr*wnoa*sin(atan(abs(hrr))))
              hrr = tan(asin(hrorder/(2.0*hrdgr*wnoa)))
              print '(" Old, new hrr, hrorder = ",2f8.3,f9.3)', &
                hrr0,hrr,hrorder
              hrrb = 10.0-200./wnoa
              hrorderb = nint(2.0*hrdgr*wnoa*sin(atan(hrrb)))
              if(hrorder/=hrorderb) then
                print '("*** Warning: hrr =",f7.3, &
                " is closer to blaze")', hrrb
                if(iwin1>0) then
                  print'("Are you sure you want hrr =",f7.3,"? (n) "$)', hrr
                  read '(a8)', yn
                  if(index(yn,'y')==0) then
                    print'("Use hrr =",f7.3,"? (y) "$)', hrrb
                    read '(a8)', yn
                    if(index(yn,'n')>0) then
                      print '(" Enter hrr "$)'
                      read *, hrr
                    else
                      hrr = hrrb
                    endif
                  endif
                endif
              endif
              if(checkall.and.(modecard/=MNONE)) wno0 = -abs(wno0)
            endif
            if(hrr0<0.) hrr = hrr0
          elseif(index(pname,'xdfl')>0) then
            read(unit=pval,fmt=*,iostat=ierr) xdfl0
          elseif(index(pname,'xdr')>0) then
            print '(" Attempt to set xdr; will be set from angle")'
          elseif(index(pname,'xdg')>0) then
            read(unit=pval,fmt=*,iostat=ierr) xdg
          elseif(index(pname,'xddgr')>0) then
            print '(" Attempt to set xddgr; will be set from mode")'
          elseif(index(pname,'mrdgr')>0) then
            read(unit=pval,fmt=*,iostat=ierr) xdmrdgr
          elseif(index(pname,'lrdgr')>0) then
            read(unit=pval,fmt=*,iostat=ierr) xdlrdgr
          elseif(index(pname,'mr0')>0) then
            read(unit=pval,fmt=*,iostat=ierr) xdmr0
          elseif(index(pname,'lr0')>0) then
            read(unit=pval,fmt=*,iostat=ierr) xdlr0
          elseif(index(pname,'x0brl')>0) then
            read(unit=pval,fmt=*,iostat=ierr) x0brl
          elseif(index(pname,'y0brl')>0) then
            read(unit=pval,fmt=*,iostat=ierr) y0brl
          elseif(index(pname,'brl')>0) then
            read(unit=pval,fmt=*,iostat=ierr) brl
          elseif(index(pname,'fftt')>0) then
            ffttort = logstr(pval,ierr)
          elseif(index(pname,'kfix')>0) then
            kfix = logstr(pval,ierr)
          elseif(index(pname,'verbose')>0) then
            verbose = logstr(pval,ierr)
            if(verbose.and.(iwin1<1)) then
              print '(" Can''t use verbose if not plotting")'
              verbose = .false.
            endif
          elseif(index(pname,'ask')>0) then
            ask = logstr(pval,ierr)
          elseif(index(pname,'pause')>0) then
            pause = logstr(pval,ierr)
            if(pause.and.(iwin1<1)) then
              print '(" Can''t use pause if not plotting")'
              pause = .false.
            endif
          elseif(index(pname,'test')>0) then
            testrun = logstr(pval,ierr)
          elseif(index(pname,'sinc')>0) then
            sincwt = logstr(pval,ierr)
          elseif(index(pname,'shiftint')>0) then
            read(unit=pval,fmt=*,iostat=ierr) (intshift(i),i=1,6)
            if(ierr/=0) then
              intshift(5) = 0
              intshift(6) = 0
              read(unit=pval,fmt=*,iostat=ierr) (intshift(i),i=1,4)
              if(ierr/=0) then
                intshift(3) = 0
                intshift(4) = 0
                read(unit=pval,fmt=*,iostat=ierr) (intshift(i),i=1,2)
              endif
            endif
          elseif(index(pname,'shift')>0) then
            doshift = logstr(pval,ierr)
          elseif(index(pname,'addwt')>0) then
            doaddwt = logstr(pval,ierr)
          elseif(index(pname,'extwt')>0) then
            if(logstr(pval,ierr)) then
              modeext = MCONWT
            else
              modeext = MUNWT
            endif
          elseif(index(pname,'weight')>0) then
            doaddwt = logstr(pval,ierr)
            if(doaddwt) then
              modeext = MCONWT
              intext(1) = 0
              intext(2) = 0
              intext(3) = 0
              intext(4) = 0
              intshift(3) = 0
              intshift(4) = 0
              intshift(5) = 0
              intshift(6) = 0
            else
              modeext = MUNWT
              intext(1) = 0
              intext(2) = 0
              intext(3) = 0
              intext(4) = 0
              intshift(3) = 0
              intshift(4) = 0
              intshift(5) = 0
              intshift(6) = 0
            endif
          elseif(index(pname,'extmode')>0) then
            if(index(pval,'none')>0) then
              modeext = MNONE
            elseif(index(pval,'unw')>0) then
              modeext = MUNWT
            elseif(index(pval,'nod')>0) then
              modeext = MNODWT
            elseif(index(pval,'conint')>0) then
              modeext = MCONINT
            elseif(index(pval,'con')>0) then
              modeext = MCONWT
            elseif(index(pval,'int')>0) then
              modeext = MINTWT
            elseif(index(pval,'avg')>0) then
              modeext = MAVGWT
            else
              print '("*** Unknown extraction mode")'
              modeext = MUNKNOWN
            endif
            if(index(pval,'int')==0) then
              intext(1) = 0
              intext(2) = 0
              intext(3) = 0
              intext(4) = 0
              intshift(3) = 0
              intshift(4) = 0
              intshift(5) = 0
              intshift(6) = 0
            endif
          elseif(index(pname,'extint')>0) then
            read(unit=pval,fmt=*,iostat=ierr) (intext(i),i=1,4)
            if(ierr/=0) then
              read(unit=pval,fmt=*,iostat=ierr) (intext(i),i=1,2)
              intext(3) = 0
              intext(4) = 0
            endif
          elseif(index(pname,'subsky')>0) then
            dosubsky = logstr(pval,ierr)
            intsky(1) = 0
            intsky(2) = 0
            intsky(3) = 0
            intsky(4) = 0
          elseif(index(pname,'skyint')>0) then
            read(unit=pval,fmt=*,iostat=ierr) (intsky(i),i=1,4)
            if(ierr/=0) then
              read(unit=pval,fmt=*,iostat=ierr) (intsky(i),i=1,2)
              intsky(3) = 0
              intsky(4) = 0
            endif
          elseif(index(pname,'illint')>0) then
            read(unit=pval,fmt=*,iostat=ierr) (intill(i),i=1,2)
            intill(3) = 0
            intill(4) = 0
          elseif(index(pname,'crossdisp')>0) then
            crossdisp = logstr(pval,ierr)
            hrr = 9.8
          elseif(index(pname,'abba')>0) then
            abba = logstr(pval,ierr)
          elseif(index(pname,'include')>0) then
            if(index(pval,'all')>0) then
              do iz = 1,mp
                wt(iz) = 1.
              enddo   
            else
              read(unit=pval,fmt=*,iostat=ierr) i1,i2
              if(ierr/=0) then
                read(unit=pval,fmt=*,iostat=ierr) i1
                i2 = i1
              endif
              if(ierr==0) then
                do iz = i1,i2
                  wt(iz) = 1.
                enddo   
              endif
            endif
          elseif(index(pname,'only')>0) then
            do iz = 1,mp
              wt(iz) = 0.
            enddo   
            read(unit=pval,fmt=*,iostat=ierr) i1,i2
            if(ierr/=0) then
              read(unit=pval,fmt=*,iostat=ierr) i1
              i2 = i1
            endif
            if(ierr==0) then
              do iz = i1,i2
                wt(iz) = 1.
              enddo   
              if(i2==i1) write(unit=ired,fmt='(i1)') i1
            endif
          elseif((index(pname,'excl')>0) &
            .or.(index(pname,'skip')>0)) then
            if(index(pval,'all')>0) then
              do iz = 1,mp
                wt(iz) = 0.
              enddo   
            else
              read(unit=pval,fmt=*,iostat=ierr) i1,i2
              if(ierr/=0) then
                read(unit=pval,fmt=*,iostat=ierr) i1
                i2 = i1
              endif
              if(ierr==0) then
                do iz = i1,i2
                  wt(iz) = 0.
                enddo   
              endif
            endif
          elseif(index(pname,'neg')>0) then
            if(index(pval,'all')>0) then
              print '(" Note: the use of ''negative = all'' to store", &
                " unsummed scans is obsolete",/ &
                " Use ''sumscan = false''")'
              do iz = 1,mp
                wt(iz) = -1.
              enddo   
            else
              read(unit=pval,fmt=*,iostat=ierr) i1,i2
              if(ierr/=0) then
                read(unit=pval,fmt=*,iostat=ierr) i1
                i2 = i1
              endif
              if(ierr==0) then
                do iz = i1,i2
                  wt(iz) = -1.
                enddo   
              endif
            endif
          elseif(index(pname,'vel')>0) then
            if(index(pval,'geo')>0) then
              call veltype(0)
            elseif(index(pval,'hel')>0) then
              call veltype(1)
            elseif(index(pval,'lsr')>0) then
              call veltype(2)
            endif
          elseif(index(pname,'sumspec')>0) then
            if(dofits.and.kbd) then
              print '(" Can''t make fitsum name in kbd mode")'
              return
            elseif(index(pval,'true')>0) then
              inquire(unit=iusumh,opened=sumopen)
              if(sumopen) then
                dosum = .true.
              else
                print'(" Sumspec not initialized.  Summing to tempsum")'
                last = index(pipedir,' ')-1
                sumfile = pipedir(1:last)//'tempsum'
                sumfits = 'tempsum.fits'
                call zerosum(ierr)
                dosum = .true.
              endif
            elseif(index(pval,'false')>0) then
              dosum = .false.
            else
              do ip = 1,70
                if(ip==70) then
                  print '(" Summing to tempsum")'
                  last = index(pipedir,' ')-1
                  sumfile = pipedir(1:last)//'tempsum'
                  sumfits = 'tempsum.fits'
                  call zerosum(ierr)
                  dosum = .true.
                  exit
                elseif(pval(ip:ip)/=' ') then
                  last = index(pipedir,' ')-1
                  sumfile = pipedir(1:last)//pval(ip:ip+30-last)
                  last = index(pval(ip:70),' ')-2
                  sumfits = pval(ip:ip+last)//'.newsum.fits'
                  call zerosum(ierr)
                  dosum = .true.
                  exit
                endif
              enddo
            endif
          elseif(index(pname,'red')>0) then
            do nskip = 1,16
              if(index(pval,' ')>1) exit
              ival = ival+1
              pval = line(ival:79)
            enddo
            ired = pval(1:32)
          elseif((index(pname,'sumscan')>0).or. &
              (index(pname,'scansum')>0)) then
            doscansum = logstr(pval,ierr)
          elseif(index(pname,'objtype')>0) then
            if(index(pval,'targ')>0) then
              objtype = 'targ'
            elseif(index(pval,'comp')>0) then
              objtype = 'comp'
            elseif(index(pval,'cal')>0) then
              objtype = 'cal'
            else
              print '(" Unknown objecttype",a20)', pval
              objtype = 'unkn'
            endif
          else
            ierr = 9
          endif
        elseif(index(line,'.')>0) then
          if(date==0.) then
            if(rdfits) then
              print '(" Enter date (yy.mmdd) "$)'
              read *, date
            else
              last = index(line,'.')+4
              cardfile = line(1:last)
              ierr = -1
              call readhdr(.false.,ierr)
              if(ierr==9) return
            endif
            call settorts(date)
          endif
          iname = index(line,'.')-1
          iseqno = iname+2
          fname = line(1:iname)
          seqno = line(iseqno:iseqno+3)
          if((checkall.or.(checkwno.and.(fname/=fnold))) &
            .and.(modecard/=MNONE)) &
            wno0 = -abs(wno0)
          fnold = fname
          if(testrun.and.(wno0>0.)) return
          call cook(fname,seqno,pipedir,fitsdir,ierr)
          if(ierr==9) then
            kbd = .true.
            ierr = 0
            return
          elseif(ierr==8) then
            ierr = 0
            return
          endif
          do iz = 1,mp
            wt(iz) = 1.
          enddo   
          ired = ' '
          baddata = .false.
          warning = ' '
          if(ierr/=0) print '(" Cook returned ierr =",i3)', ierr
          ierr = 0
        elseif(index(line,'verbose')>0) then
          if(iwin1<1) then
            print '(" Can''t use verbose if not plotting")'
            verbose = .false.
          else
            verbose = .not.verbose
            print '(" verbose = ",l1)', verbose
          endif
        elseif(index(line,'ask')>0) then
          ask = .not.ask
          print '(" ask = ",l1)', ask
        elseif(index(line,'pause')>0) then
          if(iwin1<1) then
            print '(" Can''t use pause if not plotting")'
            pause = .false.
          else
            pause = .not.pause
            print '(" pause = ",l1)', pause
          endif
        elseif(index(line,'test')>0) then
          testrun = .not.testrun
          print '(" testrun = ",l1)', testrun
        elseif(index(line,'PID')>0) then
          print '(" PID = ",a16)', pid
        elseif(index(line,'PI')>0) then
          print '(" PI = ",a16)', piname
        elseif(index(line,'flip')>0) then
          print '(" byte flip on input = ",l1)', doflip
        elseif(index(line,'flop')>0) then
          print '(" byte flip on output = ",l1)', doflop
        elseif(index(line,'pipedir')>0) then
          print '(" pipedir = ",a8)', pipedir
        elseif(index(line,'mon')>0) then
          print '(" mongo = ",l1)', domon
        elseif(index(line,'fits')>0) then
          print '(" fitsdir = ",a8)', fitsdir
          print '(" dofits = ",l1)', dofits
          print '(" rdfits = ",l1)', rdfits
        elseif(index(line,'date')>0) then
          print '(" date = ",f8.4)', date
        elseif((index(line,'cardmode')>0) &
          .or.(index(line,'flatmode')>0)) then
          if(modecard==MBLKSKY) then
            print '(" flatmode = blksky")'
          elseif(modecard==MBLKOBJ) then
            print '(" flatmode = blkobj")'
          elseif(modecard==MBLKSHINY) then
            print '(" flatmode = blkshiny")'
          elseif(modecard==MBLK) then
            print '(" flatmode = blk")'
          elseif(modecard==MSKY) then
            print '(" flatmode = sky")'
          elseif(modecard==MSHINY) then
            print '(" flatmode = shiny")'
          elseif(modecard==MNONE) then
            print '(" flatmode = none")'
          elseif(modecard==MOBJ) then
            print '(" flatmode = obj")'
          elseif(modecard==MBSBS) then
            print '(" flatmode = bsbs")'
          elseif(modecard<0) then
            print '(" flatmode = old")'
          else
            print '(" Unknown flatmode: ",i4)', modecard
          endif
        elseif(index(line,'tortmode')>0) then
          if(modetort==MBLK) then
            print '(" tortmode = blk")'
          elseif(modetort==MSKY) then
            print '(" tortmode = sky")'
          elseif(modetort==MSHINY) then
            print '(" tortmode = shiny")'
          elseif(modetort==MNONE) then
            print '(" tortmode = none")'
          elseif(modetort==MOBJ) then
            print '(" tortmode = obj")'
          elseif(modetort==MCELL) then
            print '(" tortmode = cell")'
          elseif(modetort==MOLD) then
            print '(" tortmode = old")'
          else
            print '(" Unknown tortmode: ",i4)', modetort
          endif
        elseif(index(line,'thr')>0) then
          print '(" thresh = ",f8.3)', thrfac
        elseif(index(line,'spike')>0) then
          print '(" spike = ",f8.2)', spikefac
        elseif((index(line,'stdfac')>0) &
          .or.(index(line,'noise')>0)) then
          print '(" stdfac = ",f8.2)', stdfac
        elseif(index(line,'sat')>0) then
          print '(" satval = ",f8.3)', satval/2.**16
        elseif(index(line,'dark')>0) then
          print '(" darkval = ",f10.2)', darkval
        elseif(index(line,'xnlin')>0) then
          print '(" xnlin = ",es10.2)', xnlin
        elseif(index(line,'ynlin')>0) then
          print '(" ynlin = ",es10.2)', ynlin
        elseif(index(line,'znlin')>0) then
          print '(" znlin = ",es10.2)', znlin
        elseif(index(line,'nlin')>0) then
          print '(" nlin = ",3es10.2)', xnlin,ynlin,znlin
        elseif(index(line,'dtcal')>0) then
          print '(" dtcal = ",f8.2)', dtcal
        elseif(index(line,'dtblk')>0) then
          print '(" dtcal = ",f8.2)', dtcal
        elseif(index(line,'dtsky')>0) then
          print '(" dtsky = ",f8.2)', dtsky
        elseif(index(line,'cloud')>0) then
          print '(" cloud = ",es10.2)', cloud
        elseif(index(line,'bounce')>0) then
          print '(" bounce = ",f8.2)', bounce
        elseif(index(line,'slitrot')>0) then
          print '(" slitrot = ",f8.3)', slitrot
        elseif(index(line,'lskrot')>0) then
          print '(" mrkrot = ",f8.3)', lskrot
          print '(" lrkrot = ",f8.3)', lskrot-0.004
        elseif(index(line,'xdkrot')>0) then
          print '(" xdkrot = ",f8.3)', xdkrot
        elseif(index(line,'krot')>0) then
          print '(" krot = ",f8.3)', krot
        elseif(index(line,'detrot')>0) then
          print '(" detrot = ",f8.3)', detrot
        elseif(index(line,'hrfl')>0) then
          print '(" hrfl = ",f8.2)', hrfl0
        elseif(index(line,'fred')>0) then
          print '(" fred0 = ",f8.2)', fred0
        elseif(index(line,'hrr')>0) then
          print '(" hrr = ",f8.3)', hrr
        elseif(index(line,'hrg')>0) then
          print '(" hrgamma = ",f8.3)', hrg
        elseif(index(line,'hrdgr')>0) then
          print '(" hrdgr = ",f8.4)', hrdgr
        elseif(index(line,'wno')>0) then
          print '(" wno0 = ",f9.3)', wno0
        elseif(index(line,'xdfl')>0) then
          print '(" xdfl = ",f8.2)', xdfl0
        elseif(index(line,'xdr')>0) then
          print '(" xdr =",f8.3," will be set from angle")', xdr
        elseif(index(line,'xdg')>0) then
          print '(" xdgamma = ",f8.3)', xdg
        elseif(index(line,'xddgr')>0) then
          print '(" xddgr = ",f9.6," will be set from mode")', xddgr
        elseif(index(line,'mrdgr')>0) then
          print '(" mrdgr = ",f9.6)', xdmrdgr
        elseif(index(line,'lrdgr')>0) then
          print '(" lrdgr = ",f9.6)', xdlrdgr
        elseif(index(line,'mr0')>0) then
          print '(" mr0 = ",f10.3)', xdmr0
        elseif(index(line,'lr0')>0) then
          print '(" lr0 = ",f10.3)', xdlr0
        elseif(index(line,'x0brl')>0) then
          print '(" x0brl = ",es10.2)', x0brl
        elseif(index(line,'y0brl')>0) then
          print '(" y0brl = ",es10.2)', y0brl
        elseif(index(line,'brl')>0) then
          print '(" brl = ",3es10.2)', brl,x0brl,y0brl
        elseif(index(line,'shift')>0) then
          print '(" shift = ",l1)', doshift 
          print '(" shiftint = ",6i6)', (intshift(i),i=1,6)
        elseif(index(line,'addwt')>0) then
          print '(" addwt = ",l1)', doaddwt 
        elseif(index(line,'extwt')>0) then
          if(modeext==MUNWT) then
            print '(" extwt = F")'
          else
            print '(" extwt = T")'
          endif
        elseif(index(line,'weight')>0) then
          print '(" addwt = ",l1)', doaddwt
          if(modeext==MUNWT) then
            print '(" extwt = F")'
          else
            print '(" extwt = T")'
          endif
        elseif(index(line,'extmode')>0) then
          if(modeext==MNONE) then
            print '(" extmode = none")'
          elseif(modeext==MUNWT) then
            print '(" extmode = unwt")'
          elseif(modeext==MNODWT) then
            print '(" extmode = nod")'
          elseif(modeext==MCONINT) then
            print '(" extmode = conint")'
          elseif(modeext==MCONWT) then
            print '(" extmode = con")'
          elseif(modeext==MINTWT) then
            print '(" extmode = int")'
          else
            print '(" Unknown extraction mode")'
          endif
        elseif(index(line,'extint')>0) then
          print '(" extint = ",4i6)', (intext(i),i=1,4)
        elseif(index(line,'subsky')>0) then
          print '(" subsky = ",l1)', dosubsky
        elseif(index(line,'skyint')>0) then
          print '(" skyint = ",6i6)', (intsky(i),i=1,4)
        elseif(index(line,'illint')>0) then
          print '(" illint = ",6i6)', (intill(i),i=1,4)
        elseif(index(line,'abba')>0) then
          print '(" abba = ",l1)', abba
        elseif(index(line,'check')>0) then
          print '(" checkwno, checkall =",2l2)', checkwno,checkall
        elseif(index(line,'kfix')>0) then
          print '(" kfix = ",l1)', kfix
        elseif(index(line,'fftt')>0) then
          print '(" ffttort = ",l1)', ffttort
        else
          print '("*** I don''t understand ",a40)', line
        endif

        if(ierr/=0) then
          print '("*** I don''t understand ",a40)', line
          write(iupl,'(" Error interpreting script")')
          kbd = .true.
          ierr = 9
        endif
        
      end


      subroutine settorts(date)

        real krot,lskrot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /tort0/ hrfl0,xdfl0,xdmrdgr,xdlrdgr,xdmr0,xdlr0,fred0, &
                xdkrot,lskrot
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf

        print '(" Setting tort params for date =",f8.4)', date
!	predicted or default numbers:
!	hrfl0 = 101.2
!	xdfl0 = 75.9
!	note: changed fred so nominal hrfl0 works
!	brl = 1.5
!	x0brl = 50.0
!	y0brl = 0.0
!	xdmr0 = 0.2  !  b/c xddgr is different in irobs
!	detrot = 0.080	! .083 in longslit
!	slitrot = 0.0
        irtf = .true.
!	date format: yy.mmdd
!        if(date>7.05) xdfl0 = 76.3	! don't know why this changed
        lskrot = 0.012
        xdkrot = -0.020
!        if(date<5.07) xdmr0 = 0.
        if(date<3.08) then
          hrfl0 = 100.8
          xdfl0 = 76.6
        endif
        if(date<(1.00)) then
!	I think this really means using the ZnSe focal reducer
          irtf = .false.
          fred0 = 1.94
          slitrot = -0.06
          brl = 0.0
          xdkrot = -0.015
          lskrot = 0.022
        elseif(date<(1.03)) then
          irtf = .false.
          fred0 = 1.94
          slitrot = 0.08
          brl = 0.0
          xdkrot = -0.015
          lskrot = 0.012
        elseif(date<(1.07)) then
          fred0 = 2.13
!          irtf = .false.
          slitrot = -0.01
          brl = 2.25
          x0brl = 20.
          y0brl = 0.0
          xdkrot = -0.010
          lskrot = 0.010
        elseif(date<(2.00)) then
          fred0 = 2.13
          slitrot = 0.06
          if(date<1.11145) slitrot = 0.033
          xdkrot = -0.00
        elseif(date<(2.09075)) then
          fred0 = 2.13
          slitrot = -0.02
          xdkrot = -0.02
        elseif(date<(2.12)) then
          fred0 = 2.13
          slitrot = -0.06
          xdkrot = -0.015
        elseif(date<(2.1212)) then
          fred0 = 2.13
          slitrot = -0.06
          xdkrot = -0.03
        elseif(date<(3.00)) then
          fred0 = 2.13
          slitrot = -0.06
          xdkrot = -0.025
        elseif(date<(3.08)) then
          fred0 = 2.12
          slitrot = -0.06
          xdkrot = -0.018
          if(date>3.0620) lskrot = -0.005
          x0brl = 50
        elseif(date<(5.00)) then
          slitrot = -0.06
          xdkrot = -0.02
          lskrot = 0.015
        elseif(date<(6.00)) then
          slitrot = -0.04
          xdkrot = -0.02
        elseif(date<(6.04)) then
          irtf = .false.
          slitrot = -0.08
          xdkrot = -0.02
        elseif(date<(6.09)) then
          irtf = .false.
          slitrot = -0.05
          xdkrot = -0.01
          lskrot = 0.022
        elseif(date<(6.11215)) then
          irtf = .false.
          mr0 = 0.3
          slitrot = -0.05
          xdkrot = -0.02
        elseif(date<(7.00)) then
          irtf = .false.
          mr0 = 0.3
          slitrot = -0.07
          xdkrot = -0.01
        elseif(date<(7.05)) then
          slitrot = -0.07
          xdkrot = -0.01
          y0brl = 20.
        elseif(date<(9.00)) then
          if((date>(7.09)).and.(date<(8.00))) irtf = .false.
          slitrot = -0.05
          xdkrot = -0.015
          lskrot = 0.024
        elseif(date<(10.00)) then
          slitrot = -0.08
          xdkrot = -0.015
          lskrot = 0.012 
        elseif(date<(11.00)) then
          slitrot = -0.05
          xdkrot = -0.015
          lskrot = -0.006
        elseif(date<(13.00)) then
          slitrot = -0.05
          xdkrot = -0.015
          lskrot = -0.001
        elseif(date<(13.11)) then
          slitrot = -0.07
          lskrot = 0.010 
        elseif(date<(14.00)) then
          irtf = .false.
          slitrot = -0.07
          lskrot = 0.020 
        elseif(date<(14.03)) then
!	parameters for feb14 IRTF with bad echelon chamber alignment
          hrg = 0.010
          slitrot = -0.09
          lskrot = 0.018
        elseif(date<(14.08)) then
          xdlr0 = 0.25
          slitrot = -0.08
        elseif(date<(14.09)) then
          irtf = .false.
          slitrot = -0.10
        elseif(date<(14.12)) then
          slitrot = -0.08
          lskrot = 0.014
        elseif(date<(15.00)) then
          slitrot = -0.06
          lskrot = 0.012
        elseif(date<(15.05)) then
          slitrot = -0.06
          lskrot = 0.008
        elseif(date<(17.03)) then
          slitrot = -0.06
          lskrot = 0.012
        elseif(date<(17.04)) then
          irtf = .false.
          slitrot = -0.06
          lskrot = 0.012
        elseif(date<(17.07)) then
          slitrot = -0.02
          lskrot = 0.002
        elseif(date<(18.00)) then
          slitrot = -0.06
          lskrot = 0.020
        elseif(date<(19.00)) then
          slitrot = -0.02
          lskrot = -0.001
        elseif(date<19.09) then
          slitrot = -0.02
          lskrot = 0.005
        elseif(date<22.07) then
          slitrot = -0.01
          xdkrot = 0.015
          lskrot = 0.00
        elseif(date<22.08) then
          slitrot = -0.01
          xdkrot = 0.015
          lskrot = -0.008
        elseif(date<23.00) then
          slitrot = -0.01
          xdkrot = -0.001
          lskrot = -0.008
        else
          slitrot = -0.01
          xdkrot = -0.001
          lskrot = -0.006
!	note: dates are in yy.mmdd format
        endif
        crossdisp = .true.
        krot = xdkrot

        return
      end


      function logstr(str,ierr)

        character(80) str
        character(4) str4
        logical logstr

        str4 = str(1:4)
        if((index(str4,'t')>0).or.(index(str4,'y')>0)) then
          logstr = .true.
        elseif((index(str4,'f')>0).or.(index(str4,'n')>0)) then
          logstr = .false.
        else
          ierr = 9
          logstr = .false.
        endif

        return
      end


      subroutine waitasec(pause,image,arr,nx,lx,ny)

        logical pause,image,ok
        logical pgcurs,pgband
        real arr(nx,ny)
        character(1) ch

        if(pause) then
          x = 0.
          y = 0.
          print '(" Click on a pixel for value" &
                " or hit RETURN to continue")'
!   20     ok = pgcurs(x,y,ch)
   20     ok = pgband(7,1,x,y,x,y,ch)
          if(ichar(ch)<=32) return
          if(image) then
            ll = max0(lx,ny)
            ix = x*ll+1.
            iy = y*ll+1.
            ix = min0(nx,max0(1,ix))
            iy = min0(ny,max0(1,iy))
            arri = abs(arr(ix,iy))
            if((arri>(0.1)).and.(arri<(1.e5))) then
              print '(2i8,f10.3)', ix,iy,arr(ix,iy)
            else
              print '(2i8,es10.2)', ix,iy,arr(ix,iy)
            endif
          else
            print '(2f12.3)', x,y
          endif
          goto 20
        elseif(image) then
          call sleep(2)
        else
          call sleep(1)
        endif

      end



      subroutine mswait(ms,x)

        nln = 38000
        do ims = 1,ms
          x = 0.
          do iln = 1,nln
            xiln = iln
            x = x+log(xiln)
          enddo
        enddo

      end


      subroutine cook(fname,seqno,pipedir,fitsdir,ierr)

!	routine to read in a raw set of texes data,
!	process data and write back to disk

        use ius
        use dims
        use modes

        real abeams(mx,my,mp),bbeams(mx,my,mp)
        real diffs(mx,my,mp)
        real avar(mx,my),bvar(mx,my),std(mx,my)
        real flat(mx,my),diff(mx,my)
        real spec(ms,mt),plot(mx,mt)
        real cs1(ms),cs2(ms),cx1(mx),cx2(mx),correl(3)

        character(24) fname
        character(4) seqno,ired
        character(8) pipedir,fitsdir
        character(32) rawfile,cardfile,lastcard,redfile,sumfile,posfile
        character(48) rawfits,cardfits,redfits,sumfits
        character(32) nomon
        character(60) object,feature,obsmode,instmode
        character(60) objtype,piname,pid,note,weather,warning,version
        character(60) comment
        character(8) yn
        real nodpa,lores,kmirror,krot,lskrot
        logical good,leftopen,doplot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical scan,ifobj,ifold
        real :: wnolast = -99.
        logical :: askpa = .true.

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
                piname,pid,note,weather,warning,version
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
        common /tort0/ hrfl0,xdfl0,xdmrdgr,xdlrdgr,xdmr0,xdlr0,fred0, &
                xdkrot,lskrot
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /filenames/rawfile,cardfile,lastcard,redfile,sumfile,posfile, &
                rawfits,cardfits,redfits,sumfits,ired
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,mc),sumtime
        common /iwin/ iwin1,iwin2

        save imm,wnolast,askpa
        nomon = ' '
        doplot = (iwin1>0)

!	make file names

        call makenames(fname,seqno,pipedir,fitsdir)
        call setgood(date,ierr)

        if((modecard==MNONE).and.((modetort==MNONE) &
          .or.(modetort==MOBJ).or.(modetort==MCELL))) goto 20
        if(modecard<0) then
          cardfile = lastcard
          modecard = -modecard
          ifold = .true.
        elseif(rdfits) then
          cardfile = 'tempcard'
          call fitstotx(cardfits,cardfile,ierr)
          if(ierr/=0) then
            print '(" Skipping file")'
            return
          endif
          ifold = .false.
        else
          ifold = .false.
        endif

!	if using sky from object file read object file first
        if((modecard==MOBJ).or.(modecard==MBLKOBJ)) then
          if(rdfits) then
            print '(" OBJ sky modes not implemented reading raw fits")'
            ierr = 2
            return
          endif
          jerr = -1
          call readhdr(.true.,jerr)
          np = -2
          nodon = .true.	! why?
          ierr = 0
          call readraw(rawfile,abeams,bbeams,np,ierr)
          close(unit=iuredh)
          if(dofits) close(unit=iufith)
          if(np<1) then
            print '(" Found",i4," pairs for sky in ",a16)', np,rawfile
            ierr = 2
            return
          endif
        endif

        if(modecard==MOBJ) goto 20
!	read and parse the card header

   10   if(.not.rdfits) then
          ifobj = .false.
          ierr = -1
          call readhdr(ifobj,ierr)
          close(unit=iuredh)
          if(ierr/=0) print '(" readhdr returned ierr =",i4)', ierr
          if(ierr<8) ierr = 0
        endif
        print '(2a16,f8.3,8x,a16)', &
          object,feature,waveno0,instmode
        write(iupl,'(2a16,f8.3,8x,a16)') &
          object,feature,waveno0,instmode
        if(ierr/=0) then
          if(cardfile/=lastcard) then
            print'(" Error reading ",a9,".hd  Use previous flat.")', &
              cardfile
            write(iupl, &
              '(" Error reading "a9,".hd  Use previous flat.")') &
              cardfile
            cardfile = lastcard
            ierr = 0
            goto 10
          else
            print '(" Error reading card header.  Skip file.")'
            write(iupl,'(" Error reading card header.  Skip file.")')
            ierr = 0
            return
          endif
        endif

!	read and process card data 
!	(makeflat calls testtort, which finds orders)

   20   if(modecard==MOBJ) then
!	had changed this so MOBJ meant read flat file for (blk-sky)/blk
!	but get flat from min(ab,bb) from object file
!	changing it back to reduce jan05 data
          iopt = 4
          call objcard(rawfile,abeams,bbeams,np,iopt)
        elseif(modecard==MSKY) then
          call readcard(ierr)
        elseif((modecard==MBLKOBJ).and.(index(object,'flat')/=1)) then
          if(ifold) then
            iopt = 3
            call objcard(rawfile,abeams,bbeams,np,iopt)
          else
            if((index(object,'marsflat')>0) &
              .or.(index(object,'venusflat')>0)) then
              call marscard(abeams,bbeams,ierr)
            else
              call readcard(ierr)
            endif
            iopt = 2
            call objcard(rawfile,abeams,bbeams,np,iopt)
          endif
        elseif((index(object,'marsflat')>0) &
          .or.(index(object,'venusflat')>0)) then
          call marscard(abeams,bbeams,ierr)
        elseif(modecard==MBSBS) then
          call readcard(ierr)
          if(ierr<=1) then
            do ic = 1,2
              do iy = 1,my
                do ix = 1,mx
                    cards(ix,iy,ic+2) = cards(ix,iy,ic)
                enddo   
              enddo   
            enddo   
            read(unit=seqno,fmt='(i4)') iseq
            write(unit=seqno,fmt='(i4)') iseq+1
            cardfile = 'flat.'//seqno
            call readcard(ierr)
            if(ierr<=1) then
              do ic = 1,2
                do iy = 1,my
                  do ix = 1,mx
                      cards(ix,iy,ic) = &
                        (cards(ix,iy,ic)+cards(ix,iy,ic+2))/2.
                  enddo   
                enddo   
              enddo   
            else
              print '(" Next card not found.  Set flatmode = blksky")'
              modecard = MBLKSKY
              ierr = 1
            endif
            write(unit=seqno,fmt='(i4)') iseq
            cardfile = 'flat.'//seqno
          endif    
        elseif(modecard/=MNONE) then
!	flats and fowlerflats
          call readcard(ierr)
          iopt = 0
        endif
        if(ierr>1) then
          if(.not.ifold) then
            print '(" Error reading ",a9,"  Use previous flat.")', &
              cardfile
            write(iupl,'(" Error reading ",a9,"  Use previous flat.")') &
              cardfile
            cardfile = lastcard
            ifold = .true.
            ierr = 0
            goto 10
          else
            print '(" Error reading ",a9,"  Skip file.")', cardfile
            write(iupl,'(" Error reading ",a9,"  Skip file.")') &
              cardfile
            ierr = 0
            return
          endif
        endif
        if(modecard/=MNONE) lastcard = cardfile

        call makeflat(flat,ierr)
        if(modecard/=MNONE) then
          if(ierr/=0) return
          if(verbose) then
            print '(" Plotting flat")'
            if(pause.and.ask) then
              call grey(mx,nx,ny,flat,-1.,-1.,iwin1)
            else
              call grey(mx,nx,ny,flat,0.,0.,iwin1)
            endif
            call waitasec(pause,.true.,flat,mx,nx,ny)
          endif
          if(modetort/=MNONE) then
!            call setillum(flat,ierr)
            if(pause.and.(.not.verbose)) then
              print '(" Plotting flattened card diff")'
              call grey(mx,nx,ny,cards,0.,0.,iwin1)
              call waitasec(pause,.true.,cards,mx,nx,ny)
            endif
            call tort(cards,4,ierr)
            call newillum(cards,ierr)
            if(pause.and.(.not.verbose)) then
              print '(" Plotting deskewed card diff")'
              call grey(mx,nx,ny,cards,0.,0.,iwin1)
              call waitasec(pause,.true.,cards,mx,nx,ny)
            endif
          endif
          close(unit=iuredh)

          if((wno0<=0.).or.(wno0/=wnolast)) then
            ierr = -1
            call extract(diffs,flat,std,spec,nomon,ierr)
            if(ierr/=0) print '(" extract returned ierr =",i4)', ierr
            if((ierr==8).or.(ierr==9)) return
            wnolast = wno0
          endif
        endif

        if(testrun) goto 900

        if(index(fname,'flat')>0) then
          ierr = 0
          call readhdr(.true.,ierr)
          call storespec(spec,redfile,redfits,ierr)
          goto 900
        endif

!	read and parse the raw data header while writing the reduced header

        if(rdfits) then
          rawfile = 'tempraw'
          call fitstotx(rawfits,rawfile,ierr)
          if(ierr/=0) then
            print '(" Skipping file")'
            close(unit=iufraw)
            return
          endif
        else
          jerr = 0
          call readhdr(.true.,jerr)
          if(jerr/=0) then
            print '(" Error reading data header.  Skip file.")'
            write(iupl,'(" Error reading data header.  Skip file.")')
            ierr = 0
            close(unit=iuredh)
            if(dofits) close(unit=iufith)
            return
          endif
        endif

        if(slitpa>720.) then
          print '(" slitpa not set.  Please enter slitpa (0-360): "$)'
          read *, slitpa
          slitpa = -slitpa
        endif
!	comment this out to ask for slitpa if slitpa /= nodpa
        askpa = .false.
        wavelen = 1.0e+04/waveno0
!	note: dist is now in pixels
        nodon = ((modeobs==MNOD).or.(modeobs==MMAP)) &
          .and.(((modeinst==MHIMED).and.(dist<(2.5*wavelen))) &
          .or.((modeinst==MHILOW).and.(dist<(wavelen/2.))) &
          .or.((modeinst==MMED).and.(dist<96.)) &
          .or.((modeinst==MLOW).and.(dist<96.)))
        if(nodon.and.(abs(slitpa)<360.) &
          .and.(abs(mod(abs(slitpa),180.)-mod(abs(nodpa),180.))>1.)) then
          print '("*** nodpa /= slitpa ",2f8.2)', nodpa,abs(slitpa)
          if(askpa) then
            print '(" Were we nodding along the slit? ",$)'
            read '(a8)', yn
            if(index(yn,'y')>0) then
              nodon = .true.
              print '(" Enter slitpa: ",$)'
              read *, slitpa
              slitpa = -slitpa
            else
              nodon = .false.
            endif
            askpa = .false.
          elseif(abs(mod(abs(slitpa),180.)-mod(abs(nodpa),180.))>5.) then
            nodon = .false.
            print '(" setting nodon = F")'
          endif
        endif
        scan = (modeobs==MSCAN).or.(modeobs==MMAP)
        if(scan) then
          im = 0
          modeold = modeext
          modeext = MNONE
          if((.not.dosum).or.(sumtime==0.)) then
            ierr = 0
            call zerosum(ierr)
            imm = 0
          endif
        else
          imm = 0
        endif
        np = -1

!	read raw data

  100   ierr = 0
        if(scan) then
          im = im+1
          if(((.not.doscansum).or.(wt(im)<=0.)).and.(im<=99)) then
!	reread data header to write out reduced header.im
            write(unit=ired,fmt='(".",i2.2)') im
            call makenames(fname,seqno,pipedir,fitsdir)
            if(im==1) then
              close(unit=iuredh)
              if(dofits) close(unit=iufith)
            endif
            jerr = 0
            if(.not.rdfits) call readhdr(.true.,jerr)
            if(jerr/=0) then
              print '(" Why?")'
              ierr = 0
            endif
          endif
        endif
        call readraw(rawfile,abeams,bbeams,np,ierr)
        if(scan) then
          if(wt(im)==0.) then
            ierr = -1
            goto 90
          endif
        elseif((np>=16).and.(.not.dosubsky)) then
          wt(1) = 0.
        elseif((index(obsmode,'fowler')>=1) &
          .and.(date>(7.10)).and.(date<(7.11))) then
          wt(1) = 0.
        endif
        if((ierr/=0).or.(np<=0)) then
          if(scan.and.(im>0)) then
            wt(im) = 0.
            nscan = im-1
            np = nnod
            print '(i3," complete scans found")', nscan
            goto 90
          endif
          print '(" Error reading ",a16,i4,"  Skip file.")', &
            rawfile,ierr
          write(iupl,'(" Error reading ",a16,"  Skip file.")') &
            rawfile
          ierr = 0
          close(unit=iuredh)
          if(dofits) close(unit=iufith)
          return
        endif
!        if(baddata) then
!          if(scan) goto 100
!          ierr = 0
!          close(unit=iuredh)
!          if(dofits) close(unit=iufith)
!          return
!        endif

!	despike, flat-field, calibrate, fix bads, and correct distortions

        print '(" Despiking data")'
        varfac = 1.
        if((modeobs/=MFLAT).and.(np>4)) &
          call despike(abeams,avar,varfac,np,ierr)
        if(.not.nodon) varfac = 2.
        if((modeobs/=MSTARE).and.(modeobs/=MFLAT) &
          .and.(modeobs/=MSCAN).and.(np>4)) &
          call despike(bbeams,bvar,varfac,np,ierr)
        if((modetort==MOBJ).or.(modetort==MCELL)) then
          call testtort(abeams,ierr)
          call setillum(flat,ierr)
        endif
        if((modeinst/=MCAMERA).and.(modeobs/=MSTARE) &
          .and.(modeobs/=MFLAT)) then
          if(bounce/=0.) then
            call debounce(abeams,bbeams,flat,np,ierr)
            if(scan.and.(ierr>0).and.(wt(im)>0.)) then
              if(doplot) then
                print '(" Keep scan? (y) ",$)'
                read '(a8)', yn
                if(index(yn,'n')>0) wt(im) = 0.
              else
                print '("*** Keeping scan with bad bounce")'
              endif
            endif
          endif
          if(dosubsky.and.(.not.nodon)) then
            if(.not.scan) then
              if(intsky(1)<0) then
                call cirrus(abeams,bbeams,flat,np,ierr)
              elseif(intsky(1)>0) then
                call subcorr(abeams,bbeams,flat,np,ierr)
              endif
            elseif((modeobs==MMAP).and.(wt(im)/=0.)) then
              call subcorr(abeams,bbeams,flat,np,ierr)
            endif
          endif
        endif
        ierr = 0
        if((modeobs==MSTARE)) ierr = 1
        if(modeobs==MSCAN) ierr = 2
        call diffarr(abeams,bbeams,diffs,avar,bvar,std,np,ierr)
        if((modeobs==MSCAN).and.(np<mp)) np = np+1
          if(verbose.and.(modeobs/=MFLAT)) then
            print '(" Plotting abeam(pair1)")'
            if(pause.and.ask) then
              call grey(mx,nx,ny,abeams,-1.,-1.,iwin1)
            else
              call grey(mx,nx,ny,abeams,0.,0.,iwin1)
            endif
            call waitasec(pause,.true.,abeams,mx,nx,ny)
            print '(" Plotting raw diffs(pair1)")'
            if(pause.and.ask) then
              call grey(mx,nx,ny,diffs,-1.,-1.,iwin1)
            else
              call grey(mx,nx,ny,diffs,0.,0.,iwin1)
            endif
            call waitasec(pause,.true.,diffs,mx,nx,ny)
          endif
!        if((modeobs==MFLAT).and.verbose) then
!          print '(" Plotting scatter plot of blk vs. shiny")'
!          scl = 0.
!          call scatter(nx*ny,diffs(1,1,1),diffs(1,1,3),scl,iwin1)
!          call waitasec(pause,.false.,diffs,mx,nx,ny)
!        endif
        call calibrate(diffs,flat,np,ierr)
        print '(" Cleaning diffs")'
        call clean(diffs,std,np,ierr)
!	if(verbose) then
!	  scl = 0.
!	  do ip = 1,np
!	    call grey(mx,nx,ny,diffs(1,1,ip),scl,scl,iwin1)
!	    if(pause) call waitasec(.true.,.true.,diffs(1,1,ip),mx,nx,ny)
!     	  enddo   
!	endif
        if((modeobs==MFLAT).and.verbose) then
          print '(" Plotting four flattened flat frames")'
          scl = 0.
          if(pause.and.ask) scl = -1.
          do iz = 1,4
            call grey(mx,nx,ny,diffs(1,1,iz),scl,scl,iwin1)
            call waitasec(pause,.true.,diffs(1,1,iz),mx,nx,ny)
          enddo   
        endif
        if((modeobs==MSTARE).and.(modetort==MNONE).and.verbose) then
          print '(" Plotting diffs")'
          scl = 0.
          do ip = 1,np
            call grey(mx,nx,ny,diffs(1,1,ip),scl,scl,iwin1)
            if(pause) call waitasec(ask,.true.,diffs(1,1,ip),mx,nx,ny)
          enddo   
        endif
        if(modecard/=MNONE) call calibrate(std,flat,1,ierr)
        if(modetort/=MNONE) then
          print '(" Torting diffs")'
          if(verbose) print '(" Plotting deskewed abeam")'
          call tort(abeams,1,ierr)
          if(verbose) print '(" Plotting deskewed diffs")'
          call tort(diffs,np,ierr)
          if(verbose) print '(" Plotting deskewed flattened std")'
          call tort(std,1,ierr)
          if((modeinst==MHIMED).or.(modeinst==MHILOW)) then
            slitlen = 0.5*spacing
          elseif((modeinst==MMED).or.(modeinst==MLOW)) then
            slitlen = 96.
          else
            slitlen = 0.
          endif
          nodon = (dist<slitlen)
        endif
        ierr = 0
!        if(scan.and.(im>1)) ierr = 1
        call checkarr(std,'     std',ierr)
        print '(" Mean of std =",es9.2)', amean(std)
!	obsolete.  used stdwt in common for crossdisp
!        if(stdwt) then
!          armin = arms(diffs(1,1,1))
!          do ip = 1,np
!            armsi = arms(diffs(1,1,ip))
!            if(armsi<armin) armin = armsi
!          enddo   
!          armin = 1.05*armin
!          do ip = 1,np
!            armsi = arms(diffs(1,1,ip))
!            if(armsi>armin) then
!              print '(" Skipping noisy pair ",i2,2es10.2)', &
!                ip,armsi,armin
!              wt(ip) = 0.
!            endif
!           enddo   
!        endif

   90   if(scan) then
            if(dosubsky) then
!	      if((modeobs==MSCAN).and.(intsky(1)/=0)) then
              if(modeobs==MSCAN) then
                call subscansky(diffs,std,nnod,ierr)
              elseif(nodon) then
                call submean(diffs,std,nnod,ierr)
              endif
            endif
            if((.not.doscansum).or.(wt(im)<0.)) then
!	write out an individual scan
              if((dofits).and.(wt(im)==0.)) then
                comment = 'Scan deleted from sum'
                call fithcomm('COMMENT ',comment,iufith)
                comment = 'Data flagged as bad.  Probably should ignore'
                call fithlog('BADDATA ',.true.,comment,iufith)
              endif
              ierr = -im
              if(imm==1) ierr = -1
              call storescan(diffs,flat,std,nnod, &
                redfile,redfits,beamtime,1,ierr)
            endif
!	moved this if down so bad scans are written but not summed
          if(wt(im)/=0.) then
            imm = imm+1
            call sumscan(diffs,scansum,std,np,imm,sumtime,ierr)
          endif
!	loop back to read next scan
          if(im<nscan) goto 100
!	all scans have been read from file
          print '(i3," good scans found")', imm
          if(modetort/=MNONE) then
            if(verbose) print '(" Plotting flat")'
            call tort(flat,1,ierr)
!	call extract to make wno, atmo, and noise lines
            ierr = 0
            call extract(diffs,flat,std,spec,nomon,ierr)
            if((ierr>0).and.(.not.doscansum)) then
              print '("*** Warning: ",a16, &
              " already stored.  Rerun to get wavenos right")', redfile
              print '(" Hit RETURN to continue",$)'
              read '(a8)', yn
            endif
          endif
          if(imm>0) then
            ierr = 0
            if(dosum) then
!	store scan sum over multiple files
!	this is done for every file, overwriting sumfile
              ierr = 1
              call storescan(scansum,flat,std,nnod, &
                sumfile,sumfits,sumtime,imm,ierr)
              print '(" Note: not storing individual scan files")'
            else
              if((.not.doscansum).or.(wt(im)<0.)) then
!	reread data header to write out reduced header
                ired = ' '
                call makenames(fname,seqno,pipedir,fitsdir)
                jerr = 0
                if(.not.rdfits) call readhdr(.true.,jerr)
              endif
!	store sum of scans in file (not done if summing files)
              ierr = 0
              call storescan(scansum,flat,std,nnod, &
                redfile,redfits,sumtime,imm,ierr)
            endif
          endif
          modeext = modeold
!	skip remainder of cook for scans
          goto 900
        endif

        if(modecard/=MNONE) then
          if(modetort/=MNONE) then
            if(verbose) print '(" Plotting deskewed flat")'
            call tort(flat,1,ierr)
          endif
          ierr = 0
          call checkarr(flat,'    flat',ierr)
        endif

!	correlate & shift, weight & add

        if(modetort==MNONE) then
          if(doshift) print '(" Warning: Can''t shift untorted data")'
          if(doaddwt) print '(" Warning: Can''t weight untorted data")'
          ierr = 0
!          if(modeobs==MFLAT) ierr = 1
          call adddiffs(diffs,diff,std,np,ierr)
          addtime = beamtime*np
          iwterr = 0
        else
!	  if(dosubsky.and.nodon) call submean(diffs,flat,np,ierr)
          if(dosubsky.and.nodon) call submed(diffs,flat,np,ierr)
          if(doshift) call shift(diffs,flat,std,np,ierr)
          if((modeobs==MFLAT).or.(modecard==MNONE)) then
!            ierr = 1
            call adddiffs(diffs,diff,std,np,ierr)
          else
            call wtadd(diffs,diff,flat,std,np,ierr)
          endif
          if(ierr==1) return
          iwterr = ierr
        endif

        if(iwin1>0) print '(" Plotting coadded diff")'
        if(pause.and.ask) then
          call grey(mx,nx,ny,diff,-1.,-1.,iwin1)
        else
          call grey(mx,nx,ny,diff,0.,0.,iwin1)
        endif
        call waitasec(pause,.true.,diff,mx,nx,ny)
        rmsdiff = arms(diff)
        print '(" Mean, rms of diff =",2es9.2)', amean(diff),rmsdiff
        if(verbose) then
          nyp = ny/32
          do iy = 1,nyp
            jy1 = 32*iy-31
            jy2 = 32*iy
            do ix = 1,nx
              diffsum = 0.
              do jy = jy1,jy2
                diffsum = diffsum &
                  +max(-5.*rmsdiff,min(5.*rmsdiff,diff(ix,jy)))
              enddo
              plot(ix,iy) = diffsum/(jy2-jy1+1)
            enddo
          enddo
          print '(" Plotting coadded diff collapsed by 32")'
          scale = -1.
          call perspec(mx,nyp,plot,scale,iwin2)
          call waitasec(pause,.false.,plot,mx,nx,nyp)
        endif
        if((modeobs==MSTARE).and.(modecard==MNONE)) then
          eperadu = 100.
          if(temp==(273.16)) print '(" Warning: temp not read")'
          pnut = pnu(waveno0,temp,ierr)
          aomega = pixelwd**2*(3.14159/4.)/36.
          if((modeinst==MHIMED).or.(modeinst==MHILOW)) then
            dwno = waveno0*slitval(slit,date)/(2.*abs(hrr)*hrfl0)
            rqe = (amean(diff)*eperadu)/(pnut*aomega*dwno)
            print '(" Mean RQE for black = ",es9.2, &
            " (illuminated assuming",0p,f6.0," e/ADU)")', rqe,eperadu
          elseif((modeinst==MMED).or.(modeinst==MLOW)) then
            dwno = waveno0*slitval(slit,date)/(2.*xdr*xdfl0)
            rqe = (amean(diff)*eperadu)/(pnut*aomega*dwno)
            print '(" Mean RQE for black = ",es9.2, &
            " (illuminated assuming",0p,f6.0," e/ADU)")', rqe,eperadu
          endif
        endif

!	store processed array, extract spectrum, and store spectrum

        if((modeinst==MCAMERA).or.(modeext==MNONE) &
          .or.(modetort==MNONE)) then
          if((modetort==MNONE).and.(modeext/=MNONE)) &
            print '(" Warning: Can''t extract from untorted data")'
          call storearr(diff,std,redfile,redfits,ierr)
          if(dosum.and.(iwterr==0).and.(.not.baddata)) &
            call sumarr(diff,ierr)
        elseif(crossdisp) then
          ierr = 0
          call extract(diff,flat,std,spec,redfile,ierr)
          call storespec(spec,redfile,redfits,ierr)
          if(dosum.and.(iwterr==0).and.(.not.baddata)) &
            call sumspec(spec,ierr)
        else
          ierr = 0
          call extract(diff,flat,std,spec,redfile,ierr)
          call storearr(diff,std,redfile,redfits,ierr)
          if(dosum.and.(iwterr==0).and.(.not.baddata)) &
            call sumarr(diff,ierr)
        endif

  900   if(doplot.and.(modetort/=MNONE).and.(modecard/=MNONE)) then
          if((modeinst==MHIMED).or.(modeinst==MHILOW)) then
            do ix = 1,nx
              sumc = 0.
              sumc2 = 0.
              sumi = 0.
              do iy = ny/4+1,3*ny/4
                cardi = cards(ix,iy,1)
                if(cardi/=0.) then
                  sumc = sumc+cardi
                  sumc2 = sumc2+cardi**2
                  sumi = sumi+1.
                endif
              enddo
              sumi = max(1.,sumi)
              sumc = sumc/sumi
              stddev = sqrt(abs(sumc2/sumi-sumc**2))
              cmx = sumc+5.*stddev
              cmn = sumc-5.*stddev
              sumc = 0.
              sumi = 0.
              do iy = 1,ny
                cardi = cards(ix,iy,1)
                if((cardi/=0.).and.(cardi>cmn).and.(cardi<cmx)) then
                  sumc = sumc+cardi
                  sumi = sumi+1.
                endif
              enddo
              sumc = sumc/max(1.,sumi)
              do iy = 1,ny
                cardi = cards(ix,iy,1)
                if(cardi/=0.) then
                  if(cardi<cmn) then
                    cards(ix,iy,1) = cmn-sumc
                  elseif(cardi>cmx) then
                    cards(ix,iy,1) = cmx-sumc
                  else
                    cards(ix,iy,1) = cards(ix,iy,1)-sumc
                  endif
                endif
              enddo
            enddo
          else
            do iy = 1,ny
              sumc = 0.
              sumi = 0.
              do ix = 1,nx
                if(cards(ix,iy,1)/=0.) then
                  sumc = sumc+cards(ix,iy,1)
                  sumi = sumi+1.
                endif
              enddo
              sumc = sumc/max(1.,sumi)
              do ix = 1,nx
                if(cards(ix,iy,1)/=0.) &
                  cards(ix,iy,1) = cards(ix,iy,1)-sumc
              enddo
            enddo
            iy1 = 1
            iy4 = ny
            do iy = 1,ny
              if(illx(iy)>0) then
                if(iy1==1) then
                  iy1 = iy
                else
                  iy4 = iy
                endif
              endif
            enddo
            iy2 = (iy1+iy4-1)/2
            iy3 = iy4-(iy2-iy1)
            do ix = 1,nx
              sumc = 0.
              do iy = iy1,iy2
                sumc = sumc+cards(ix,iy,1)
              enddo
              cx1(ix) = sumc
              sumc = 0.
              do iy = iy3,iy4
                sumc = sumc+cards(ix,iy,1)
              enddo
              cx2(ix) = sumc
            enddo
            do ishft = -1,1
              sumc = 0.
              do ix = 9,nx-8
                sumc = sumc+cx1(ix)*cx2(ix+ishft)
              enddo
              correl(ishft+2) = sumc
            enddo
            dcdsh = (correl(3)-correl(1))/2.
            ddcdsh = correl(1)-2.*correl(2)+correl(3)
            if(ddcdsh<abs(dcdsh)) then
              dsh = sign(1.,dcdsh)
            else
              dsh = -dcdsh/ddcdsh
            endif
            angle = -dsh/(iy3-iy1)
            print '(" krot,angle,recommend: ",3f8.4)', &
              krot,angle,krot-angle/2.
          endif
          if(testrun) then
            print '(" Plotting sky/blk-mean")'
            call grey(mx,nx,ny,cards,0.,0.,iwin1)
            call waitasec(.true.,.true.,cards,mx,nx,ny)
          endif
        endif

        if(rdfits)then
          print '(" Finished ",a24)', rawfits
        else
          print '(" Finished ",a16)', rawfile
        endif
        inquire(unit=iuredh,opened=leftopen)
        if(leftopen) then
!          print '(" Closing reduced header ",a32)', redfile
          close(unit=iuredh)
          if(dofits) close(unit=iufith)
        endif
        if(rdfits) close(unit=iufraw)

        return
      end


      subroutine makenames(fname,seqno,pipedir,fitsdir)

        use ius

        character(24) fname
        character(4) seqno,ired
        character(8) pipedir,fitsdir
        character(24) fitsdate
        character(32) rawfile,cardfile,lastcard,redfile,sumfile,posfile
        character(48) rawfits,cardfits,redfits,sumfits
        character(60) line
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        real nodpa,lores,kmirror,krot

        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /filenames/rawfile,cardfile,lastcard,redfile,sumfile,posfile, &
                rawfits,cardfits,redfits,sumfits,ired
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl

        lenn = index(fname,' ')-1
        if(lenn<0) lenn = 24
        rawfile = fname(:lenn)//'.'//seqno
        lenr = index(rawfile,' ')-1
        cardfile = 'flat.'//seqno
        lenc = index(cardfile,' ')-1
        lenp = index(pipedir,' ')-1
        ilen = index(ired,' ')-1
        if((ilen<0).or.(ilen>4)) ilen = 4
        if(ilen>0) then
          if((lenp+lenn)>25) lenn = 25-last-ilen
          redfile = &
          pipedir(:lenp)//fname(:lenn)//'.'//seqno//ired(:ilen)
        else
          if((lenp+lenn)>25) lenn = 25-lenp
          redfile = pipedir(:lenp)//fname(:lenn)//'.'//seqno
        endif
        iyr = int(date)
        iday = nint(10000.*(date-iyr))
        if(iday<200) then
          iyr = iyr-1
          write(unit=fitsdate,fmt='("TX",i2.2,"B",i4.4,".")') iyr,iday
        elseif(iday<800) then
          write(unit=fitsdate,fmt='("TX",i2.2,"A",i4.4,".")') iyr,iday
        else
          write(unit=fitsdate,fmt='("TX",i2.2,"B",i4.4,".")') iyr,iday
        endif
        lend = index(fitsdate,' ')-1
        if(lend<0) lend = 24
        lenf = index(fitsdir,' ')-1
        if(lenf<0) lenf = 8
!	I replaced rawfile(:lenr) and cardfile(:lenf) with seqno
        if(rdfits) then
!	in this case these files are used for input
          rawfits = fitsdate(:lend)//seqno//'.raw.fits'
          cardfits = fitsdate(:lend)//seqno//'.flt.fits'
          if(lenf>0) then
            fitsdate = fitsdir(:lenf)//fitsdate(:lend)
            lend = min(lenf+lend,24)
          endif
          if(ilen>0) then
            redfits = &
              fitsdate(:lend)//seqno//ired(:ilen)//'.red.fits'
          else
            redfits = fitsdate(:lend)//seqno//'.red.fits'
          endif
          if(index(sumfits,'newsum')>0) then
            lens = index(sumfits,'newsum')-2
            sumfits = fitsdate(:lend)//sumfits(:lens)//'.red.fits'
          endif
          rawfile = 'tempraw'
          cardfile = 'tempcard'
        elseif(dofits) then
          if(lenf>0) then
            fitsdate = fitsdir(:lenf)//fitsdate(:lend)
            lend = min(lenf+lend,24)
          endif
          rawfits = fitsdate(:lend)//seqno//'.raw.fits'
          cardfits = fitsdate(:lend)//seqno//'.flt.fits'
          if(ilen>0) then
            redfits = &
              fitsdate(:lend)//seqno//ired(:ilen)//'.red.fits'
          else
            redfits = fitsdate(:lend)//seqno//'.red.fits'
          endif
          if(index(sumfits,'newsum')>0) then
            read(unit=seqno,fmt=*) minseq
            maxseq = minseq
            nline = 0
            do
              read(unit=iups,fmt='(a60)',iostat=ierr) line
              if((ierr>0).or.(index(line,'end')>0)) then
                print '(" Error reading pipescript to get summed files")'
                exit
              endif
              nline = nline+1
              if(index(line,'storesum')>0) exit
              do ichar = 1,8
                if(line(1:1)==' ') then
                  line(1:59) = line(2:60)
                else
                  exit
                endif
              enddo
              ipound = index(line,'#')
              if(ipound==1) then
                cycle
              elseif(ipound==0) then
                ipound = 20
              else
                ipound = ipound-1
              endif
              if(index(line(1:ipound),'=')>0) cycle
              iseq = index(line(1:ipound),'.')+1
              if(iseq<3) cycle
              read(unit=line(iseq:ipound),fmt='(i4)',iostat=ierr) newseq
              if(ierr/=0) cycle
              minseq = min(minseq,newseq)
              maxseq = max(maxseq,newseq)
            enddo
            write(unit=sumfits,fmt='(i4.4,"-",i4.4)') minseq,maxseq
            sumfits = fitsdate(:lend)//sumfits(1:9)//'.red.fits'
            do iline = 1,nline
              backspace(unit=iups)
            enddo
            ierr = 0
          endif
        else
          rawfits = 'raw.fits'
          cardfits = 'flt.fits'
          redfits = 'red.fits'
          sumfits = 'sum.fits'
        endif

        return
      end


      subroutine setgood(date,ierr)
!	initialize good array and flag known bad pixels

        use dims

        logical good
        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)

        do iy = 1,ny
          do ix = 1,nx
            good(ix,iy) = .true.
          enddo   
        enddo   

        do iy = 133,134
          do ix = 23,29
            good(ix,iy) = .false.
          enddo   
        enddo   
        do iy = 132,135
          do ix = 23,26
            good(ix,iy) = .false.
          enddo   
        enddo   
        good(185,136) = .false.

        ierr = 0
        return
      end



      subroutine setone(arr,nz,ierr)
!	set an arr to 1

        use dims

        real arr(mx,my,nz)
        common /nn/ nx,ny,nc,norder,ns,nt

        do iz = 1,nz
          do iy = 1,ny
            do ix = 1,nx
              arr(ix,iy,iz) = 1.
            enddo   
          enddo   
        enddo   

        ierr = 0
        return
      end


      function amean(arr)
!	return mean of illuminated pixels in arr

        use dims

        real arr(mx,my)
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)

        sum = 0.
        isum = 0
        do iy = 1,my
          do ix = 1,mx
            if(illum(ix,iy)>0) then
              if((0.*arr(ix,iy))==0.) then
                sum = sum+arr(ix,iy)
                isum = isum+1
              endif
            endif
          enddo   
        enddo   
        if(isum>0) then
          amean = sum/isum
        else
!	 print '(" illum = 0 for all pixels")'
          amean = 0.
        endif

        return
      end


      function arms(arr)
!	return rms of illuminated pixels in arr

        use dims
        real arr(mx,my)
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)

        sum = 0.
        isum = 0
        do iy = 1,my
          do ix = 1,mx
            if(illum(ix,iy)>0) then
              if((0.*arr(ix,iy))==0.) then
                sum = sum+arr(ix,iy)**2
                isum = isum+1
              endif
            endif
          enddo   
        enddo   
        if(isum>0) then
          arms = sqrt(sum/isum)
        else
!	 print '(" illum = 0 for all pixels")'
          arms = 0.
        endif

        return
      end


      subroutine readhdr(ifobj,ierr)
!	read and interpret a raw header, and copy to red header and fits.hd

        use ius
        use dims
        use modes
        use consts
        character(26), parameter :: capchar = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        character(26), parameter :: lowchar = 'abcdefghijklmnopqrstuvwxyz'

        character(32) rawfile,cardfile,lastcard,redfile,sumfile,posfile
        character(48) rawfits,cardfits,redfits,sumfits
        character(4) ired
        character(36) rawhd,cardhd,redhd,sumhd
        character(80) line,val
        character(60) object,feature,obsmode,instmode
        character(60) objtype,piname,pid,note,weather,warning,version
        character(60) comment,NOC,posline
        logical ifobj,leftopen
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical doplot
        real nodpa,lores,kmirror,krot,lskrot
        real :: pa = 0.,oldpa = 0.
        character(8) yn
        logical hdcopy

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
                piname,pid,note,weather,warning,version
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
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /tort0/ hrfl0,xdfl0,xdmrdgr,xdlrdgr,xdmr0,xdlr0,fred0, &
                xdkrot,lskrot
        common /filenames/rawfile,cardfile,lastcard,redfile,sumfile,posfile, &
                rawfits,cardfits,redfits,sumfits,ired
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime
        common /iwin/ iwin1,iwin2

        save oldpa

        doplot = (iwin1>0)
        NOC = ' '
        ra = 15.*12+34/4.+56.78/240.
        equinox = 2000.0
        hdcopy = (ierr==0).and.(.not.rdfits)

!	open raw, reduced, and fits header files

        if(ifobj) then

          inquire(unit=iurawh,opened=leftopen)
          if(leftopen) then
            print '("*** Raw header was left open")'
            close(unit=iurawh)
          endif
          namelen = index(rawfile,' ')-1
          if(namelen<0) namelen = 32
          rawhd = rawfile(:namelen)//'.hd'
          open(unit=iurawh,file=rawhd,status='old',iostat=ierr)
          if(ierr/=0) then
            print '("*** Error: ",a24," does not exist.")', rawhd
            write(iupl,'("!!! Error: ",a24," does not exist.")') rawhd
            print '(" Skipping file. ", &
              " Hit RETURN to continue or q to quit ",$)'
            read '(a8)', yn
            if(index(yn,'q')>0) call quit
            ierr = 9
            return
          endif

          inquire(unit=iuredh,opened=leftopen)
          if(leftopen) close(unit=iuredh)
          if(hdcopy.and.(.not.rdfits)) then
            namelen = index(redfile,' ')-1
            if(namelen<0) namelen = 32
            redhd = redfile(:namelen)//'.hd'
            open(unit=iuredh,file=redhd)
          else
            open(unit=iuredh,status='scratch')
          endif

        else

          inquire(unit=iurawh,opened=leftopen)
          if(leftopen) then
            print '("*** Raw flat header was left open")'
            close(unit=iurawh)
          endif
          namelen = index(cardfile,' ')-1
          if(namelen<0) namelen = 32
          cardhd = cardfile(:namelen)//'.hd'
          open(unit=iurawh,file=cardhd,status='old',iostat=ierr)
          if(ierr/=0) then
            print '("*** Error: ",a24," does not exist.")', cardhd
            write(iupl,'("!!! Error: ",a24," does not exist.")') cardhd
            print '(" Skipping file.  Hit RETURN to continue",$)'
            read '(a8)', yn
            ierr = 9
            return
          endif
          inquire(unit=iuredh,opened=leftopen)
          if(leftopen) close(unit=iuredh)
          open(unit=iuredh,status='scratch')

        endif

        if(dofits) then
          inquire(unit=iufith,opened=leftopen)
!	fits.hd might have been left open after writing raw.fits
          if(leftopen) close(unit=iufith)
          open(unit=iufith,file='fits.hd')
        endif

!	read, write, and parse lines

        ierr = 0
        gain = 1.0
        obstime = 1.0
        lorder = 0
  100   read(unit=iurawh,fmt='(a79)',err=190,end=200) line
        do ich = 1,8
          icap = index(capchar,line(ich:ich))
          if(icap>0) line(ich:ich) = lowchar(icap:icap)
        enddo
        if(ifobj) write(iuredh,'(a79)') line
        if(line(1:8)=='object  ') then
          read(unit=line,fmt='(10x,a16)',err=190,end=190) object
          if(dofits) call fithchar('OBJECT  ',object,NOC,iufith)
          if((index(piname,'none')<1).and.(index(piname,' ')>1)) then
             if(ifobj) write(iuredh,'("PI      = ",a16)') piname
             if(dofits) call fithchar('PI      ',piname,NOC,iufith)
          endif
!	parameters set in pipescript
          if(index(pid,' ')>1) then
            if(ifobj) write(iuredh,'("PID     = ",a16)') pid
            comment = 'IRTF program ID'
            if(dofits) call fithchar('PID     ',pid,comment,iufith)
          endif
          if(index(weather,' ')>1) then
            if(ifobj) write(iuredh,'("weather = ",a60)') weather
            if(dofits) call fithchar('WEATHER ',weather,NOC,iufith)
          endif
          if(index(warning,' ')>1) then
            if(ifobj) write(iuredh,'("warning = ",a60)') warning
            if(dofits) call fithchar('WARNING ',warning,NOC,iufith)
          endif
          if(ifobj) then 
            if(index(objtype,' ')>1) then
              write(iuredh,'("objtype = ",a16)') objtype
              if(index(objtype,'targ')>0) then
                comment = 'primary science target'
              elseif(index(objtype,'comp')>0) then
                comment = 'telluric comparison source'
              elseif(index(objtype,'cal')>0) then
                comment = 'calibration source'
              else
                comment = 'unknown object type'
              endif
              if(dofits) call fithchar('OBJTYPE ',objtype,comment,iufith)
            endif
          endif
        elseif(line(1:8)=='feature ') then
          read(unit=line,fmt='(10x,a16)',err=190,end=190) feature
          comment = 'spectral feature'
          if(dofits) call fithchar('FEATURE ',feature,comment,iufith)
        elseif(line(1:8)=='obsmode ') then
          read(unit=line,fmt='(10x,a16)',err=190,end=190) obsmode
          comment = NOC
          if(index(obsmode,'flat')>0) then
            modeobs = MFLAT
            comment = 'blackbody flat-field source'
          elseif(index(obsmode,'nod')>0) then
            modeobs = MNOD
            comment = 'AB nod observation'
          elseif(index(obsmode,'fowler')>0) then
            if(dist/=0.) then
              modeobs = MNOD
              line = "obsmode = nod "
              comment = 'nod mode with Fowler sampling'
            elseif(dosubsky) then
              modeobs = MSCAN
              line = "obsmode = scan "
              comment = 'zero scan step with Fowler sampling'
            else
              modeobs = MSTARE
              line = "obsmode = stare "
              comment = 'zero nod throw with Fowler sampling'
            endif
            if(ifobj) write(iuredh,'(a79)') line
          elseif(index(obsmode,'chop-nod')>0) then
            modeobs = MCHOPNOD
          elseif(index(obsmode,'chop')>0) then
            modeobs = MCHOP
          elseif(index(obsmode,'stare')>0) then
            modeobs = MSTARE
          elseif(index(obsmode,'dark')>0) then
            modeobs = MSTARE
          elseif(index(obsmode,'map')>0) then
            modeobs = MMAP
          elseif(index(obsmode,'fscan')>0) then
            modeobs = MSCAN
          elseif(index(obsmode,'scan')>0) then
            modeobs = MSCAN
          else
            modeobs = MUNKNOWN
            print '("^G*** Unknown obs mode")'
            print '(" Hit RETURN to continue",$)'
            read '(a8)', yn
            write(iupl,'("!!! Unknown obs mode")')
          endif
          if(dosum.and.(modeobs>=MSCAN).and.(.not.doscansum)) then
            print '(" Summing scans over files without summing within files", &
              " is not allowed")'
            print '(" First reduce with sumscan = sumspec = true")'
            print '(" Then reduce with sumscan = sumspec = false")'
            print '(" Set sumscan = true? (y) "$)'
            read '(a8)', yn
            if(index(yn,'n')>0) stop
            doscansum = .true.
          endif
          if(dofits) call fithchar('OBSMODE ',obsmode,comment,iufith)
          if(.not.irtf) then
            if(modeobs==MFLAT) then
              comment = 'flat'
            elseif((modeobs==MUNKNOWN).or.(modeobs==MSTARE)) then
              comment = 'eng'
            else
              comment = 'object'
            endif
            if(dofits) call fithchar('OBSTYPE ',comment,NOC,iufith)
            if(modeobs==MFLAT) then
              comment = 'partnerCal'
            elseif((modeobs==MUNKNOWN).or.(modeobs==MSTARE)) then
              comment = 'eng'
            else
              comment = 'science'
            endif
            if(dofits) call fithchar('OBSCLASS',comment,NOC,iufith)
          endif
          if((index(obsmode,'fscan')>0).and.(index(posfile,' ')>1)) then
            open(unit=16,file=posfile)
            do
              read(unit=16,fmt='(a16)',iostat=ierr) posline
              if(ierr/=0) exit
              write(iuredh,'("scanstep= ",a16)') posline
              comment = 'fscan step arcsec E N'
              if(dofits) call fithchar('SCANSTEP',posline,comment,iufith)
            enddo
            close(unit=16)
          endif
        elseif(line(1:8)=='instmode') then
          if(verbose) print '(a40)', line(1:40)
          read(unit=line,fmt='(10x,a16)',err=190,end=190) instmode
          if(index(instmode,'hi-med')>0) then
            modeinst = MHIMED
            comment = 'echelon x echelle'
          elseif(index(instmode,'med')>0) then
            modeinst = MMED
            comment = 'echelle only'
          elseif(index(instmode,'hi-lo')>0) then
            modeinst = MHILOW
            comment = 'echelon x lo-res'
          elseif(index(instmode,'lo')>0) then
            modeinst = MLOW
            comment = 'lo-res only'
          elseif(index(instmode,'cam')>0) then
            modeinst = MCAMERA
            comment = 'lo-res face-on'
          elseif(index(instmode,'pup')>0) then
            modeinst = MCAMERA
            comment = 'pupil viewing camera'
          else
            modeinst = MUNKNOWN
            print '("*** Unknown inst mode")'
            print '(" Hit RETURN to continue",$)'
            read '(a8)', yn
            write(iupl,'("!!! Unknown inst mode")')
            comment = 'unknown instrument mode'
          endif
          if((modeinst==MMED).or.(modeinst==MHIMED)) then
            xddgr = xdmrdgr
          elseif((modeinst==MLOW).or.(modeinst==MHILOW)) then
            xddgr = xdlrdgr
          endif
          if(((modeinst==MHILOW).or.(modeinst==MHIMED)) &
            .and.(.not.crossdisp)) then
            crossdisp = .true.
            hrr = 9.8
!            lskrot = krot
            if(xdkrot/=0.) then
              krot = xdkrot
            else
              print '("***Why hasn''t xdkrot been set?")'
              xdkrot = -.020
              krot = -.020
            endif
            detrot = .080
            print '(" Setting crossdisp,krot = ",l1,f7.3)', crossdisp,krot
          elseif(((modeinst==MLOW).or.(modeinst==MMED)) &
            .and.crossdisp) then
            crossdisp = .false.
            hrr = 9.8
            xdkrot = krot
            if(lskrot/=0.) then
              krot = lskrot
            else
              print '("***Why hasn''t lskrot been set?")'
!	in case we really want lskrot = 0
!              lskrot = .012
              krot = lskrot
            endif
            detrot = .083
            print '(" Setting crossdisp,krot = ",l1,f7.3)', crossdisp,krot
          endif
          if(modeinst==MMED) then
            krot = lskrot
          elseif(modeinst==MLOW) then
            krot = lskrot-.004
          endif
          if(dofits) call fithchar('INSTMODE',instmode,comment,iufith)
        elseif(line(1:8)=='romtable') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) irom
          rdrst = (irom==9)
          fowler = (irom==0)
          rdhalf = (irom==12)
          if(rdrst) then
            if(darkval<(2.**15)) darkval = 15.*(2.**12)
            if(xnlin/=0.) then
              print '("*** non-linearity correction not implemented", &
                " in read-reset mode")'
              xnlin = 0.
            endif
          elseif(darkval>(2.**15)) then
            darkval = 0.
          endif
          comment = 'specifies readout mode'
          if(dofits) call fithint('ROMTABLE',irom,comment,iufith)
        elseif(line(1:8)=='waveno0 ') then
          if(verbose) print '(a40)', line(1:40)
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) waveno0
          if(waveno0<=350.) waveno0 = 1000.
          comment = 'planned central wavenumber'
          if(dofits) call fithreal('WAVENO0 ',waveno0,comment,iufith)
        elseif(line(1:8)=='order   ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) xdorder
          lorder = nint(xdorder)
          comment = 'cross-dispersion grating order number'
          if(dofits) call fithint('ORDER   ',lorder,comment,iufith)
        elseif(line(1:8)=='temp    ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) temper
          temper = temper+dtcal
          comment = 'ambient temperature (K)'
          if((temper>(-15.)).and.(temper<35.)) temp = temper+273.16
          if((temp<250.).or.(temp>310.)) then
            print '(" Implausible temperature.")'
            write(iupl,'(" Implausible temperature.")')
            comment = 'substituted for implausible value'
            temp = 273.16
          else
            temp = temper
          endif
          if(dofits) call fithreal('TEMPER  ',temp,comment,iufith)
        elseif((line(1:8)=='temp_det').or.(line(1:8)=='dewartem')) then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) tempdet
          comment = 'detector temperature (unknown units)'
          if(dofits) call fithreal('TEMP_DET',tempdet,comment,iufith)
        elseif(line(1:8)=='humidity') then
          if(index(line,'%')>0) then
            last = index(line,'%')-1
            read(unit=line(11:last),fmt=*,err=191,end=191) humid
          else
            read(unit=line,fmt='(10x,f12.0)',err=191,end=191) humid
          endif
          comment = 'dome level humidity; not related to pwv'
          if(dofits) call fithreal('HUMIDITY',humid,comment,iufith)
        elseif(line(1:8)=='nod     ') then
          iche = index(line,'E')
          ichw = max(iche,index(line,'W'))
          ichn = index(line,'N')
          if((iche<12).or.(ichn<iche+2)) goto 190
          val = line(10:iche-1)
          read(unit=val,fmt=*,err=190,end=190) dra
          val = line(ichw+1:ichn-1)
          read(unit=val,fmt=*,err=190,end=190) ddec
          dist = sqrt(dra**2+ddec**2)
          if(dofits) call fithreal('NOD     ',dist,line(11:70),iufith)
          nodpa = DEGRAD*atan2(dra,ddec)
          if(nodpa<0.) nodpa = nodpa+360.
        elseif((line(1:8)=='inttime ') &
                .or.(line(1:8)=='frametim')) then
          if(verbose) print '(a40)', line(1:40)
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) frtime
          if(frtime<=0.) frtime = 1.0
          comment = 'integration time / frame'
          if(dofits) call fithreal('FRAMETIM',frtime,comment,iufith)
        elseif(line(1:8)=='obstime ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) obstime
          if(obstime<=0.) obstime = 1.0
          comment = 'calculated observation time (s)'
          if(dofits) call fithreal('OBSTIME ',obstime,comment,iufith)
        elseif(line(1:8)=='tottime ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) tottime
          if(tottime<=0.) tottime = 1.0
          comment = 'calculated total clock time (s)'
          if(dofits) call fithreal('TOTTIME ',tottime,comment,iufith)
        elseif(line(1:8)=='echelle ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) echelle
          comment = 'echelle grating angle (deg)'
          if(dofits) call fithreal('ECHELLE ',echelle,comment,iufith)
          echelle = echelle-xdmr0
          if((modeinst==MMED).or.(modeinst==MHIMED)) &
            xdr = tan(echelle/DEGRAD)
        elseif(line(1:8)=='lores   ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) lores
          comment = 'lo-res grating angle (deg)'
          if(dofits) call fithreal('LORES   ',lores,comment,iufith)
          lores = lores-xdlr0
          if((modeinst==MLOW).or.(modeinst==MHILOW)) &
            xdr = tan(lores/DEGRAD)
        elseif(line(1:8)=='kmirror ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) kmirror
          comment = 'K-mirror angle'
          if(dofits) call fithreal('KMIRROR ',kmirror,comment,iufith)
        elseif(line(1:8)=='slit    ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) slit
          comment = 'slit wheel angle'
          if(dofits) call fithreal('SLIT    ',slit,comment,iufith)
        elseif(line(1:8)=='filter  ') then
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) filter
          comment = 'filter wheel angle'
          if(dofits) call fithreal('FILTER  ',filter,comment,iufith)
        elseif(line(1:8)=='gain    ') then
          if(verbose) print '(a40)', line(1:40)
          read(unit=line(11:20),fmt=*,err=190,end=190) gain
          if((date>(6.11)).and.(date<(6.1112))) gain = gain/2.
          if(gain==0.) gain = 1.0
          comment = 'electronic gain'
          if(dofits) call fithreal('GAIN    ',gain,comment,iufith)
        elseif(line(1:8)=='bias    ') then
          if(verbose) print '(a40)', line(1:40)
          read(unit=line(11:20),fmt=*,err=190,end=190) bias
          comment = 'detector bias voltage (V)'
          if(dofits) call fithreal('BIAS    ',bias,comment,iufith)
        elseif(line(1:8)=='pressure') then
          read(unit=line(11:20),fmt=*,err=190,end=190) pressure
          comment = 'LN2 pumping line pressure (Torr)'
          if(dofits) call fithreal('PRESSURE',pressure,comment,iufith)
        elseif(line(1:8)=='pixelwd ') then
!	no longer included in header
          read(unit=line,fmt='(10x,f12.0)',err=190,end=190) pixelwd
          comment = 'pixel width (um?)'
          if(dofits) call fithreal('PIXELWD ',pixelwd,comment,iufith)
        elseif(line(1:8)=='nframe  ') then
          if(verbose) print '(a40)', line(1:40)
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nframe
          comment = 'number of frames coadded in hardware'
          if(dofits) call fithint('NFRAME  ',nframe,comment,iufith)
        elseif(line(1:8)=='nsum    ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nsum
          comment = 'number of frames coadded in software'
          if(dofits) call fithint('NSUM    ',nsum,comment,iufith)
        elseif(line(1:8)=='nwrite  ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nwrite
          comment = 'number of frames written per nod phase'
          if(dofits) call fithint('NWRITE  ',nwrite,comment,iufith)
        elseif(line(1:8)=='nscan   ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nscan
          comment = 'number of scans'
          if(dofits) call fithint('NSCAN   ',nscan,comment,iufith)
        elseif(line(1:8)=='nnod    ') then
          if((modeobs==MSCAN).or.(modeobs==MMAP)) then
            if(index(obsmode,'fowler')>0) then
              read(unit=line,fmt='(10x,i8)',err=190,end=190) nnod
            else
              read(unit=line,fmt='(10x,i8)',err=190,end=190) nscan
              comment = 'number of scans'
              if(dofits) call fithint('NSCAN   ',nscan,comment,iufith)
            endif
          else
            read(unit=line,fmt='(10x,i8)',err=190,end=190) nnod
            comment = 'number of nod pairs'
            if(dofits) call fithint('NNOD    ',nnod,comment,iufith)
            if((modeobs==MSCAN).or.(modeobs==MMAP)) then
              if(index(obsmode,'fowler')==0) nscan = nnod
            elseif(nnod>mp) then
              print '(" nnod =",i3," > mp =",i3, &
                "  Reading only through mp")', nnod,mp
              write(iupl,'(" nnod =",i3," > mp =",i3, &
                "  Reading only through mp")') nnod,mp
              nnod = mp
              write(unit=comment,fmt='("pipe could read only",i4)') mp
              if(dofits) call fithcomm('COMMENT ',comment,iufith)
            endif
          endif
        elseif((line(1:8)=='nsteps  ').or.(line(1:8)=='npoints ')) &
          then
          if(((modeobs==MSCAN).or.(modeobs==MMAP)) &
            .and.(index(obsmode,'fowler')==0)) then
            read(unit=line,fmt='(10x,f10.0)',err=190,end=190) points
            nnod = points
            comment = 'number of points in scan'
            if(dofits) call fithint('NPOINTS ',nnod,comment,iufith)
          endif
        elseif(line(1:8)=='nsky    ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nsky
          comment = 'number of extra points on sky'
          if(dofits) call fithint('NSKY    ',nsky,comment,iufith)
          if(((modeobs==MSCAN).or.(modeobs==MMAP)) &
            .and.(index(obsmode,'fowler')==0)) then
            if((intsky(2)>nnod+nsky).or.(intsky(4)>nnod+nsky)) print &
              '("***WARNING: skyint is set outside of scan length")'
          endif
        elseif(line(1:8)=='nspec   ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nspec
          ny = nspec
          comment = 'array size in spectral dimension'
          if(dofits) call fithint('NSPEC   ',nspec,comment,iufith)
        elseif(line(1:8)=='nspat   ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nspat
          nx = nspat
          comment = 'array size in spatial dimension'
          if(dofits) call fithint('NSPAT   ',nspat,comment,iufith)
        elseif(line(1:8)=='bitpix  ') then
          read(unit=line,fmt='(10x,i8)',err=190,end=190) nbitpix
        elseif(line(1:8)=='time    ') then
          read(unit=line,fmt='(10x,i2,1x,i2,1x,f3.0)', &
            err=191,end=191) ih,im,s
          ut = ih+im/60.+s/3600.
          comment = 'UTC'
          if(dofits) call fithchar('TIME-OBS',line(11:70),comment,iufith)
        elseif(line(1:8)=='date-obs') then
          read(unit=line,fmt='(10x,i2,1x,i2,1x,i2)', &
            err=109,end=191) iyr,imo,idy
          line = 'date-obs= 20'//line(11:68)
          goto 110
  109     read(unit=line,fmt='(10x,i4,1x,i2,1x,i2)', &
            err=191,end=191) iyr,imo,idy
          iyr = iyr-2000
  110     if(date==0.) then
            date = iyr+(imo+idy/100.)/100.
            if((date>30.).or.(date<0.)) then
              print '(" date bad ",f12.4,a30)', date,val
              stop
            endif
          endif
          if(dofits) then
            comment = 'yyyy-mm-dd observation date'
            if(irtf) then
              call fithchar('DATE-OBS',line(11:70),comment,iufith)
            else
              call fithchar('DATE    ',line(11:70),comment,iufith)
            endif
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
            write(unit=val,fmt='(i4,"-",i2.2,"-",i2.2)') jyr,jmo,jdy
            comment = 'yyyy-mm-dd public release date'
            call fithchar('DATE-REL',val,comment,iufith)
          endif
        elseif(line(1:8)=='ra      ') then
          read(unit=line,fmt='(10x,i2,1x,i2,1x,f5.2)', &
            err=191,end=191) ih,im,s
          ra = 15.*ih+im/4.+s/240.
          if(dofits) call fithreal('RA      ',ra,line(11:21),iufith)
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
          if(dofits) call fithreal('DEC     ',dec,line(11:21),iufith)
        elseif(line(1:8)=='targra  ') then
          read(unit=line,fmt='(10x,i2,1x,i2,1x,f5.2)', &
            err=191,end=191) ih,im,s
          tra = 15.*ih+im/4.+s/240.
          if(dofits) call fithchar('TARG_RA ',line(11:70),NOC,iufith)
        elseif(line(1:8)=='targdec ') then
          if(line(11:11)=='-') then
            read(unit=line,fmt='(10x,i3,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            tdec = id-im/60.-s/3600.
          elseif(line(11:11)=='+') then
            read(unit=line,fmt='(11x,i2,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            tdec = id+im/60.+s/3600.
          else
            read(unit=line,fmt='(10x,i2,1x,i2,1x,f5.2)', &
              err=191,end=191) id,im,s
            tdec = id+im/60.+s/3600.
          endif
          if(dofits) call fithchar('TARG_DEC',line(11:70),NOC,iufith)
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
          if(dofits) call fithreal('AZIMUTH ',az,line(11:21),iufith)
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
          if(dofits) call fithreal('ELEVATIO',el,line(11:21),iufith)
        elseif(line(1:8)=='equinox ') then
          read(unit=line,fmt='(10x,f7.2)',err=191,end=191) equinox
          if(equinox==0.) then
            equinox = 2000.+iyr+dayn(iyr,imo,idy,0.)/365.
          elseif((equinox<1900.).or.(equinox>2050.)) then
            equinox = 2000.
          endif
        elseif(line(1:8)=='ha      ') then
          if(index(line,':')>0) then
            if(line(11:11)=='-') then
              read(unit=line,fmt='(10x,i3,1x,i2,1x,f5.2)', &
                err=193,end=191) id,im,s
              ha = id-im/60.-s/3600.
            elseif(line(11:11)=='+') then
              read(unit=line,fmt='(11x,i2,1x,i2,1x,f5.2)', &
                err=193,end=191) id,im,s
              ha = id+im/60.+s/3600.
            else
              read(unit=line,fmt='(10x,i2,1x,i2,1x,f5.2)', &
                err=193,end=191) id,im,s
              ha = id+im/60.+s/3600.
            endif
          else
            read(unit=line,fmt='(10x,f7.4)',err=191,end=191) ha
          endif
          comment = 'degrees'
          if(dofits) call fithreal('HA      ',15.0*ha,comment,iufith)
        elseif(line(1:8)=='airmass ') then
          read(unit=line,fmt='(10x,f7.3)',err=191,end=191) airmass
          comment = 'at start of integration'
          if(dofits) call fithreal('AIRMASS ',airmass,comment,iufith)
        elseif(line(1:8)=='instrpa ') then
          read(unit=line,fmt='(10x,f7.3)',err=191,end=191) pa
          comment = 'instrument position angle'
          if(dofits) call fithreal('INSTRPA ',pa,comment,iufith)
        elseif(line(1:8)=='instraa ') then
          read(unit=line,fmt='(10x,f7.3)',err=191,end=191) aa
          if(dofits) call fithreal('INSTRAA ',aa,NOC,iufith)
        elseif(line(1:8)=='telescop') then
          if(index(line,'McD')>0) then
            irtf = .false.
            efl0 = 12.*270.
          elseif(index(line,'IRTF')>0) then
            irtf = .true.
            efl0 = 12.7*300.
          elseif(index(line,'Gemini')>0) then
            irtf = .false.
            if(date<(6.04)) then
              efl0 = 10.6*800.
            elseif(date<(6.08)) then
              efl0 = 10.7*800.
            else
              efl0 = 11.5*800.
            endif
          elseif((date<(6.0)).or.(date>(8.0))) then
            irtf = .true.
            print '(" Unknown telescope")'
            efl0 = 12.*300.
            slitpa = 0.
          else
            irtf = .false.
            print '(" Unknown telescope")'
            efl0 = 12.*800.
            slitpa = 0.
          endif
          if(dofits) call fithchar('TELESCOP',line(11:70),NOC,iufith)
          comment = 'slit position angle (deg)'
!	I think pa is instrpa from the header
          if(.not.irtf) then
!	on Gemini slitpa = instrpa+180
            slitpa = mod(pa+180.,360.)
            if(dofits) call fithreal('SLITPA  ',slitpa,comment,iufith)
          elseif(slitpa<0.) then
!	slitpa has been read from pipescript
            if(dofits) call fithreal('SLITPA  ',-slitpa,comment,iufith)
          elseif(slitpa>720.) then
!	slitpa is still set to the default
            continue
          elseif(pa/=oldpa) then
!	instrpa has changed
            slitpa = slitpa-oldpa+pa
            slitpa = mod(slitpa,360.)
            if(slitpa<0.) slitpa = slitpa+360.
            print '(" Changing slitpa to",f8.3)', slitpa
            if(dofits) call fithreal('SLITPA  ',slitpa,comment,iufith)
          else
!	slitpa was calculated from nod
            if(dofits) call fithreal('SLITPA  ',slitpa,comment,iufith)
          endif
          oldpa = pa
          if(ifobj.and.(slitpa<720.)) &
            write(iuredh,'("slitpa  = ",f12.4)') abs(slitpa)
        elseif(line(1:7)=='instrum') then
          comment = 'TEXES'
          if(dofits) call fithchar('INSTRUME',comment,NOC,iufith)
        elseif(line(1:8)=='gemprgid') then
!	where did this come from?
!          dofits = (index(line,'ENG')<=0).and.(redfits/='red.fits')
          if(dofits) call fithchar(line(1:8),line(11:70),NOC,iufith)
        else
          if(dofits) call fithchar(line(1:8),line(11:70),NOC,iufith)
        endif
        goto 100

  192   if(line(11:11)=='-') then
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
        if(dofits) call fithreal('AZIMUTH ',az,line(11:21),iufith)
        goto 100
  193   if(line(11:11)=='-') then
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
        if(dofits) call fithreal('ELEVATIO',el,line(11:21),iufith)
        goto 100

  190   ierr = ierr+1
  191   if(ifobj) then
          print '(" Error reading line from ",a16/a80)', rawhd,line
          write(iupl,'(" Error reading line from ",a16/a80)') rawhd,line
        else
          print '(" Error reading line from ",a16/a80)', cardhd,line
          write(iupl,'(" Error reading line from ",a16/a80)')cardhd,line
        endif
        goto 100

  200   close(unit=iurawh)

        if(((ra-tra)>(0.5)).and.((ra-tra)<(0.64)) &
          .and.(equinox==(1950.))) then
          print '(" Changing equinox for ra, dec to 2000.")'
          equinox = 2000.
          comment = 'equinox for RA, DEC.  TARG equinox is 1950.'
        else
          comment = NOC
        endif
        if(dofits) call fithreal('EQUINOX ',equinox,comment,iufith)

        day = dayn(iyr,imo,idy,ut)
        dummy = 15.*12+34/4.+56.78/240.
        if(ifobj) then
          if(abs(ra-dummy)>(0.01)) then
            if(equinox/=(2000.0)) call precess(ra,dec,equinox,2000.0)
            call vels(day,ra,dec)
          else
            print '(" Dummy coordinates in header.", &
              "  Can''t calculate v_Earth")'
          endif
        endif

!	calculate wavelength dependence of fred and efl
!	note: fred includes effect of shifting the xd parab (if any)
!	note: efl0 is set elsewhere by date and telescope
!       xnkbr = sqrt(2.3613-3.115e+04/waveno0**2-5.86e+08/waveno0**4)
        wavel2 = (1.0e+04/waveno0)**2
        xnkbr = sqrt(2.36197+0.02293/wavel2 &
          -0.1767/((60.61e-04*waveno0)**2-1.0) &
          -2.0622/((87.72e-04*waveno0)**2-1.0))
        xncdte = sqrt(1.0+6.19779/(1.0-0.10053/wavel2) &
          -3.2244/(5279.5/wavel2-1.0))
        xnznse = sqrt(1.0+4.45814/(1.0-0.04035/wavel2) &
          +0.46722/(1.0-0.15317/wavel2)-2.89566/(2221.82/wavel2-1.0)) &
          -0.013
!	index at 10um
        xnkbr0 = 1.52714
        xncdte0 = sqrt(1.0+6.19779/(1.0-0.10053e-02)-3.2244/51.795)
        xnznse0 = sqrt(1.0+4.45814/(1.0-0.04035e-02) &
          +0.46722/(1.0-0.15317e-02)-2.89566/21.2182)-0.013
        if(date<(1.03)) then
!	ZnSe focal reducer
!          fred = 1.93*((waveno0-320.)/(waveno0-300.))/(680./700.)
!	adjust fred to get fls right
!          fred = 1.96*((waveno0-320.)/(waveno0-300.))/(680./700.)
!          fred0 = 1.96
!	now the old number seems better
          if(fred0>2.0) fred0 = 1.93
          powrat = (xnznse-1.0)/(xnznse0-1.0)
          fmag = 1.0-powrat*(1.0-1.0/fred0)
          fred = 1.0/fmag
          if(.not.ifobj) nc = 3
        else
!	KBr focal reducer
!	not sure what old formula assumed.
!	try a new one assuming lens-det dist is constant
!	adjust fred to get fls right
          if(fred0<2.0) fred0 = 2.12
          fmag0 = 1.0/fred0
          powrat = (xnkbr-1.)/(xnkbr0-1.)
          fmag = 1.0-powrat*(1.0-fmag0)
          fred = 1.0/fmag
!         pow = .525*(xnkbr-1.)/(xnkbr0-1.)
!         fred = 1./(1.-pow)
          if(.not.ifobj) nc = nnod
        endif
        if(irtf) then
          flfor0 = 90.
          flfore = flfor0*(xnkbr0-1.0)/(xnkbr-1.0)
          flcdte0 = 270.
          flcdte = flcdte0*(xncdte0-1.0)/(xncdte-1.0)
          d1 = 120.	! slit to KBr lens
          d12 = 360.	! distance between lenses
          fp = 7642.	! primary mirror fl
          fs = 665.44	! secondary mirror -fl
          dps = 7022.22	! nominal mirror separation
          ftel = 1.0/(1.0/fp-1.0/fs+dps/(fp*fs))  ! nominal tel fl
!	I think my old formulas assumed d12 is fixed and d01 varies
!	but we refocus the telescope so it is d12 that varies
          fl1 = flfor0
          fl2 = flcdte0
!	try a formula based on central rays entering dewar parallel to axis
          xmfor0 = (d12/fl2-1.0)*(1.0-d1/fl1+d1/(d12-fl2))
!	calculate the change in focus position and its effect on IRTF pltscl
          df2 = d12-1.0/(1.0/fl1-1.0/d1)
          dfoc = df2/(df2/fl2-1.0)
          ddps = dfoc/214.
!	fractional change in telescope fl
          dftel0 = ftel*ddps/(fp*fs)
          fl1 = flfore
          fl2 = flcdte
          xmfore = (d12/fl2-1.0)*(1.0-d1/fl1+d1/(d12-fl2))
          df2 = d12-1.0/(1.0/fl1-1.0/d1)
          dfoc = df2/(df2/fl2-1.0)
          ddps = dfoc/214.
          dftel = ftel*ddps/(fp*fs)
          forerat = (xmfore/xmfor0)*(1.0+dftel0-dftel)
!	this might be right without the CdTe lens
!	  xmfore = (dfore-flfore)/flfore
!	  xmfor0 = (dfore-flfor0)/flfor0
!	this looks interesting.  it gives same answer as above
!         xmfore = flfore*flcdte/((d01-flfore)*(d12-flcdte)-d01*flfore)
!         xmfor0 = flfor0*flcdte/((d01-flfor0)*(d12-flcdte)-d01*flfor0)
!	empirical correction to agree with Tommy's measurements
!         forerat = (xmfore/xmfor0)**0.7
        else
!	this seems to apply before installing current IRTF foreoptics
!	and at Gemini
          forerat = 1.
        endif
!	not sure when or why this changed
        if(date>(16.0)) fred = 1.01*fred
        hrfl = hrfl0/fred
        xdfl = xdfl0/fred
  	efl = efl0*forerat/fred
        pltscl = pixelwd/(efl*4.848e-06)
        slitwid = slitval(slit,date)/(efl0*4.848e-06)
        omegap = pltscl*slitwid*(4.848e-06)**2
        if(ifobj) then
          write(iuredh,'("pltscl  = ",f12.4)') pltscl
          write(iuredh,'("slitwid = ",f12.4)') slitwid
          if(dofits) then
            comment = 'arcsec / pixel after focal reduction'
            call fithreal('PLTSCL  ',pltscl,comment,iufith)
            comment = 'deg / pixel along slit'
            call fithreal('CDELT2  ',pltscl/3600.,comment,iufith)
            comment = 'slit width (arcsec)'
            call fithreal('SLITWID ',slitwid,comment,iufith)
          endif
          dist = dist/pltscl
        endif

        if(index(obsmode,'fowler')>0) then
          nframe = nframe/2
          frtime = (52./77.)*frtime/nframe
        endif
        beamtime = frtime*nframe*nsum*nwrite

        if(ifobj.and.(modeobs==MSCAN)) then
          if(date<(1.03)) then
            nsky = 0
            write(iuredh,'("nsky    = ",i8)') nsky
            comment = 'number of extra points on sky'
            if(dofits) call fithint('NSKY    ',nsky,comment,iufith)
          elseif(date<(1.07)) then
            nsky = (nnod+9)/10
            nsky = max0(3,nsky)
            write(iuredh,'("nsky    = ",i8)') nsky
            comment = 'number of extra points on sky'
            if(dofits) call fithint('NSKY    ',nsky,comment,iufith)
          endif
          nnod = nnod+nsky
        endif
        if(nnod>mp) then
          print '(" nnod =",i3," > mp =",i3, &
                "  Arrays will overflow")', nnod,mp
          write(iupl,'(" nnod =",i3," > mp =",i3, &
                "  Arrays will overflow")') nnod,mp
          ierr = -1
          close(unit=iuredh)
          close(unit=iufith)
          return
        endif

        goto 250
  240   continue
          print '(" lores, kmirror = ",2f8.2)', lores,kmirror
          print '(" Enter instmode ",$)'
          read '(a16)', instmode
          if(index(instmode,'hi-med')>0) then
            modeinst = MHIMED
            comment = 'echelon x echelle (corrected in pipe)'
            xddgr = xdmrdgr
          elseif(index(instmode,'med')>0) then
            modeinst = MMED
            comment = 'echelle only (corrected in pipe)'
            xddgr = xdmrdgr
          elseif(index(instmode,'hi-lo')>0) then
            modeinst = MHILOW
            comment = 'echelon x lo-res (corrected in pipe)'
            xddgr = xdlrdgr
          elseif(index(instmode,'lo')>0) then
            modeinst = MLOW
            comment = 'lo-res only (corrected in pipe)'
            xddgr = xdlrdgr
          elseif(index(instmode,'cam')>0) then
            modeinst = MCAMERA
            comment = 'lo-res face-on (corrected in pipe)'
          elseif(index(instmode,'pup')>0) then
            modeinst = MCAMERA
            comment = 'pupil viewing camera (corrected in pipe)'
          else
            print '(" instmode options:  hi-med, med, hi-lo, lo, cam")'
            goto 240
          endif
          if(dofits) call fithchar('INSTMODE',instmode,comment,iufith)
  250   if((modeinst==MHIMED).or.(modeinst==MMED)) then
          iorder = nint(2.0*xddgr*sin(echelle/DEGRAD)*waveno0)
          if((iorder-lorder)/=0) then
            print '(" Warning: header and calculated order disagree:", &
              2i4/3es10.3)', lorder,iorder,xddgr,echelle,waveno0
            print '(" Using calculated order ",i3)', iorder
            lorder = iorder
          endif
          wno = lorder/(2.0*xddgr*sin(echelle/DEGRAD))
          if(wno0==0.) then
            continue
          elseif(wno0==-2.) then
            wno0 = wno
            print '(" Setting wno0 = ",f10.3)', -wno0
          elseif(wno0==-1.) then
            wno0 = waveno0
            print '(" Setting wno0 = ",f10.3)', -wno0
          elseif(wno0==1.) then
            wno0 = -waveno0
            print '(" Setting wno0 = ",f10.3)', -wno0
          elseif(wno0==2.) then
            wno0 = -wno
            print '(" Setting wno0 = ",f10.3)', -wno0
          elseif(abs(wno-abs(wno0))>wno/1000.) then
            print '(" Warning: echelle and wno0 disagree",3f10.3)', &
              echelle,wno,abs(wno0)
            write(iupl,'("!!! Warning: echelle and wno0 disagree", &
              3f10.3)') echelle,wno,abs(wno0)
            if(doplot) then
              print '(" Change wno0? (n) ",$)'
              read '(a8)', yn
              if(index(yn,'y')>0) then
                print '(" Set wno0 = waveno0 = ",f10.3,"? (y) "$)', waveno0
                read '(a8)', yn
                if(index(yn,'n')==0) then
                  wno0 = waveno0
                else
                  print '(" Enter wno0 (0 to run setwno) "$)'
                  read *, wno0
                endif
              else
                print '(" Adjust mr0? (y) ",$)'
                read '(a8)', yn
                if(index(yn,'n')==0) then
                  xdorder = iorder
                  xdang0 = DEGRAD*asin(xdorder/(2.0*xddgr*abs(wno0)))
                  xdmr00 = xdmr0
                  xdmr0 = echelle+xdmr0-xdang0
                  if((xdmr0<0.0).or.(xdmr0>0.4)) then
                    print '(" Are you sure you want mr0 = ",f8.3,"? "$)', &
                      xdmr0
                    read '(a8)', yn
                    if(index(yn,'n')>0) then
                      print '(" Enter mr0: ",$)'
                      read *, xdmr0
                      echelle = echelle+xdmr00-xdmr0
                    else
                      echelle = xdang0
                    endif
                  else
                    echelle = xdang0
                  endif
                  print '(" mr0 = ",f8.3)', xdmr0
                  write(iupl,'(" mr0 = ",f8.3)') xdmr0
                endif
              endif
            endif
          elseif((abs(wno-abs(waveno0))>wno/100.) &
            .or.(doplot.and.(abs(wno-abs(waveno0))>wno/200.))) then
            print '("*** Warning: echelle and waveno0 disagree", &
              3f10.2)', echelle,wno,waveno0
            print '(" waveno0 is probably wrong.  set waveno0 = ",f10.3)', wno
            waveno0 = wno
          endif
          if(lores>(-12.)) then
            print '("*** Warning: lores and instmode disagree")'
            write(iupl,'("!!! Warning: lores and instmode disagree")')
            print '(" Enter corrected instmode? (n) ",$)'
            read '(a8)', yn
            if(index(yn,'y')>0) goto 240
          endif
        elseif((modeinst==MHILOW).or.(modeinst==MLOW)) then
          iorder = nint(2.0*xddgr*sin(lores/DEGRAD)*waveno0)
          if(lorder*(iorder-lorder)/=0) then
            print '(" Warning: header and calculated order disagree:", &
              2i4/3es10.3)', lorder,iorder,xddgr,lores,waveno0
            lorder = iorder
          endif
          wno = lorder/(2.0*xddgr*sin(lores/DEGRAD))
          if(wno0==0.) then
            continue
          elseif(wno0==-2.) then
            wno0 = wno
            print '(" Setting wno0 = ",f10.3)', -wno0
          elseif(wno0==-1.) then
            wno0 = waveno0
            print '(" Setting wno0 = ",f10.3)', -wno0
          elseif(wno0==1.) then
            wno0 = -waveno0
            print '(" Setting wno0 = ",f10.3)', -wno0
          elseif(wno0==2.) then
            wno0 = -wno
            print '(" Setting wno0 = ",f10.3)', -wno0
          elseif((abs(wno-abs(wno0))>wno/200.) &
            .or.(doplot.and.(abs(wno-abs(wno0))>wno/500.))) then
            print '(" Warning: lores and wno0 disagree", &
              3f10.3)', lores,wno,abs(wno0)
            write(iupl,'("!!! Warning: lores and wno0 disagree", &
              3f10.3)') lores,wno,abs(wno0)
            if(doplot) then
              print '(" Change wno0? (n) ",$)'
              read '(a8)', yn
              if(index(yn,'y')>0) then
                print '(" Set wno0 = waveno0 = ",f10.3," (y) "$)', waveno0
                read '(a8)', yn
                if(index(yn,'n')==0) then
                  wno0 = waveno0
                else
                  print '(" Enter wno0 (0 to run setwno) "$)'
                  read *, wno0
                endif
              else
                print '(" Will reset lr0")'
              endif
            endif
          elseif(abs(wno-waveno0)>(waveno0/200.)) then
            print '(" Warning: lores and waveno0 disagree", &
              3f10.2)', lores,wno,waveno0
            write(iupl,'(" Warning: lores and waveno0 disagree", &
              3f10.2)') lores,wno,waveno0
            print '(" Derived order, waveno: ",i2,f10.2)', iorder,wno
          endif
          if(lores<0.) then
            print '("*** Warning: lores and instmode disagree")'
            write(iupl,'("!!! Warning: lores and instmode disagree")')
            print '(" Enter corrected instmode? (n) ",$)'
            read '(a8)', yn
            if(index(yn,'y')>0) goto 240
          endif
        elseif((modeinst==MCAMERA).and.(abs(lores)>(2.)))then
          print '("*** Warning: lores and instmode disagree")'
          write(iupl,'("!!! Warning: lores and instmode disagree")')
          print '(" Enter corrected instmode? (n) ",$)'
          read '(a8)', yn
          if(index(yn,'y')>0) goto 240
        endif
        if((modeinst==MHIMED).or.(modeinst==MHILOW)) then
          if(kmirror>(3.)) then
            print '("*** Warning: kmirror and instmode disagree")'
            write(iupl,'("!!! Warning: kmirror and instmode disagree")')
            print '(" Enter corrected instmode? (n) ",$)'
            read '(a8)', yn
            if(index(yn,'y')>0) goto 240
          endif
        elseif((modeinst==MMED).or.(modeinst==MLOW)) then
          if((kmirror<(42.)).or.(kmirror>(48.))) then
            print '("*** Warning: kmirror and instmode disagree")'
            write(iupl,'("!!! Warning: kmirror and instmode disagree")')
            print '(" Enter corrected instmode? (n) ",$)'
            read '(a8)', yn
            if(index(yn,'y')>0) goto 240
          endif
        endif
        if(modeinst/=MCAMERA) then
!	which wno is most reliable for calculating xdr?
          wnoa = abs(wno0)
!	I would like to get xdr from middle of the array
!	but I don''t yet know how many orders are on the array
          if((modeinst==MMED).or.(modeinst==MHIMED)) then
            sinang = sin(echelle/DEGRAD)
          else
            sinang = sin(lores/DEGRAD)
          endif
          wno = lorder/(2.0*xddgr*sinang)
          if((abs(wno-wnoa)/(wno+wnoa)<0.001) &
            .or.(abs(waveno0-wnoa)/(waveno0+wnoa)<0.001)) then
            sinang = lorder/(2.0*xddgr*wnoa)
          elseif(abs(waveno0-wno)/(waveno0+wno)>0.001) then
            print '(" waveno0,wno(xdang),wno0 disagree",3f9.3)', &
              waveno0,wno,wnoa
            if(doplot) then
              print '(" Enter wno for xdr calculation: ",$)'
              read *, wno
            else
              print '(" Using wno(xdang) = ",f8.3)', wno
            endif
            sinang = lorder/(2.0*xddgr*wno)
          endif
          xdr = sinang/sqrt(1.-sinang**2)
        endif

        if(ifobj.and.hdcopy) then
          write(iuredh,'("beamtime= ",f11.6)') beamtime
          write(iuredh,'("pipeline= fife.F90 version ",a8)') version
          if(dofits) then
            comment = 'frtime*nframe*nsum*nwrite'
            call fithreal('BEAMTIME',beamtime,comment,iufith)
            comment = ' fife.F90 version '//version(1:8)
            call fithchar('PIPELINE',comment,NOC,iufith)
          endif
          ierr = 0
          call writeparms(iuredh,iufith,ierr)
          if(dosum) then
            if(sumtime==0.) then
!	copy red header to sum header
              close(unit=iuredh)
              open(unit=iuredh,file=redhd)
              inquire(unit=iusumh,opened=leftopen)
              if(leftopen) then
                print '("*** Sum header was left open")'
                backspace(unit=iusumh)
                read(unit=iusumh,fmt='(a79)',iostat=ierr) line
                if(ierr==0) print '(a80)', line
                close(unit=iusumh)
              endif
              if(rdfits) then
                open(unit=iusumh,status='scratch')
              else
                namelen = index(sumfile,' ')-1
                if(namelen<0) namelen = 32
                sumhd = sumfile(:namelen)//'.hd'
                open(unit=iusumh,file=sumhd)
                do
                  read(unit=iuredh,fmt='(a79)',iostat=ierr) line
                  if(ierr/=0) exit
                  write(iusumh,'(a79)') line
                enddo
                backspace(unit=iuredh)
                ierr = 0
              endif

              if(dofits) then
!	copy fits.hd to fitsum.hd
                close(unit=iufith)
                open(unit=iufith,file='fits.hd')
                inquire(unit=iufish,opened=leftopen)
                if(leftopen) then
                  print '("*** Fitsum header was left open")'
                  read(unit=iufish,fmt='(a79)',iostat=ierr) line
                  if(ierr==0) print '(a80)', line
                  close(unit=iufish)
                endif
                open(unit=iufish,file='fitsum.hd')

                do
                  read(unit=iufith,fmt='(a78)',iostat=ierr) line
                  if(ierr/=0) exit
                  write(iufish,'(a78)') line
                enddo
                backspace(unit=iufith)
                ierr = 0
              endif

            endif
            write(iusumh,'("rawfile = ",a24)') rawfile
            comment = rawfile
            if(dofits) call fithchar('RAWFILE ',comment,NOC,iufish)
          endif
        elseif(dofits) then
          comment = 'End of flat header'
          call fithcomm('COMMENT ',comment,iufith)
        endif
        if(dofits) then
          if(index(note,' ')>0) &
            call fithcomm('COMMENT ',note,iufith)
          if(baddata) then
            comment = 'Data flagged as bad.  Probably should ignore'
            call fithlog('BADDATA ',.true.,comment,iufith)
!          else
!            comment = 'Data not flagged as bad'
!            call fithlog('BADDATA ',.false.,comment,iufith)
          endif
        endif

        if(.not.(hdcopy.and.ifobj)) then
          close(unit=iuredh)
          close(unit=iufith)
!       else leave redhd and fits.hd open
        endif

      end


      function slitval(slit,date)
!	convert slit from degrees to cm

        dimension val(24,6)
        data val/0.015,0.035,0.045,0.055,0.100,1.000, &
                 0.015,0.020,0.025,0.035,0.055,0.000, &
                 12*0.000, &
                 0.025,0.035,0.012,0.015,0.020,0.025, &
                 0.035,0.045,0.055,1.000,0.040,0.020, &
                 0.012,0.025,0.020,0.020,0.025,0.035, &
                 6*0.000, &
                 0.025,0.035,0.025,0.020,0.020,0.025, &
                 0.012,0.020,0.040,1.000,0.055,0.045, &
                 0.035,0.025,0.020,0.015,0.012,0.035, &
                 6*0.000, &
                 0.024,0.024,0.024,0.024,0.024,0.018, &
                 0.018,0.018,0.018,1.000,0.015,0.018, &
                 0.024,0.036,0.054,0.024,0.018,0.012, &
                 6*0.000, &
                 24*0.000, &
                 24*0.000/

        if(date<(1.00)) then
          if(slit<(170.)) then
            islit = nint(slit/22.5)-2
          else
            islit = nint(slit/30.)
          endif
          if((islit>=1).and.(islit<=12)) then
            slitval = val(islit,1)
          else
            print '(" ** Warning: slit angle out of range")'
            slitval = 0.01
          endif
        elseif(date<(2.00)) then
          islit = nint(slit/20.)+1
          if((islit>=1).and.(islit<=18)) then
            slitval = val(islit,2)
          else
            print '(" ** Warning: slit angle out of range")'
            slitval = 0.01
          endif
        elseif(date<(2.08)) then
!	was this slit wheel ever used?
          islit = nint(slit/20.)+1
          if((islit>=1).and.(islit<=18)) then
            slitval = val(islit,3)
          else
            print '(" ** Warning: slit angle out of range")'
            slitval = 0.01
          endif
        else
!	was this installed in sep02?
          islit = nint(slit/20.)+1
          if((islit>=1).and.(islit<=18)) then
            slitval = val(islit,4)
          else
            print '(" ** Warning: slit angle out of range")'
            slitval = 0.01
          endif
        endif

        return
      end



      function dayn(iy,im,id,ut)
!	calculate day number from J2000.0

        dimension iday(12)
        data iday/0,31,59,90,120,151,181,212,243,273,304,334/

        leap = (iy+4.+(iday(im)+id-60)/365.25)/4.
        dayn = iy*365.25+iday(im)+id+leap+ut/24.-1.5

        return
      end


      subroutine vels(day,ra,dec)
!	calculate Earth's and Sun's motion away from object

        use ius
        use dims
        use consts

        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        real krot
        character(60) comment

        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime

        save ivel

        slong = 280.454+0.985647*day
        ilong = slong/360.
        slong = slong-360.*ilong
        sanom = 357.528+0.985600*day
        ianom = sanom/360.
        sanom = (sanom-360.*ianom)/DEGRAD
        slong = (slong+1.915*sin(sanom)+0.020*sin(2.*sanom))/DEGRAD
        soblq = 23.439/DEGRAD
        sdist = 1.00014-0.01671*cos(sanom)-0.00014*cos(2.*sanom)
        gdot = 0.985600/DEGRAD
        rdot = 0.01671*sin(sanom)*gdot
        sldot = (0.985647 &
          +gdot*(1.915*cos(sanom)+0.040*cos(2.*sanom)))/DEGRAD
        xdot = rdot*cos(slong)-sldot*sin(slong)*sdist
        ydot = rdot*cos(soblq)*sin(slong) &
          +sldot*cos(soblq)*cos(slong)*sdist
        zdot = tan(soblq)*ydot

        today = 2000.0+day/365.25
        call precess(ra,dec,2000.0,today)
        ra = ra/DEGRAD
        dec = dec/DEGRAD
        xobj = cos(ra)*cos(dec)
        yobj = sin(ra)*cos(dec)
        zobj = sin(dec)

        vgeo = 1731.46*(xdot*xobj+ydot*yobj+zdot*zobj)
        vsol = 19.5*(cos(30./DEGRAD)*yobj-sin(30./DEGRAD)*zobj)
        velsr = vgeo+vsol

        print '(" Radial velocity of Earth, Sun:",2f7.2)', vgeo,vsol
        if(.not.rdfits) then
          write(iuredh,'("vehelio = ",f7.2)') vgeo
          write(iuredh,'("velsr   = ",f7.2)') velsr
        endif
        if(dofits) then
          comment = 'heliocentric motion of Earth from object(km/s)'
          call fithreal('VEHELIO ',vgeo,comment,iufith)
          comment = 'LSR motion of Earth from object (km/s)'
          call fithreal('VELSR   ',velsr,comment,iufith)
        endif

        if(ivel==1) then
          radvel = vgeo
        elseif(ivel==2) then
          radvel = velsr
        else
          radvel = 0.
        endif

        return

      entry veltype(itype)

        ivel = itype

        return

      end


      subroutine precess(ra,dec,y1,y2)
!	precess coordinates from equinox y1 to y2

        use consts

        xm(t) = 2.23617e-02*t+6.77e-06*t**2+1.75e-07*t**3
        xn(t) = 9.71717e-03*t-2.07e-06*t**2-2.02e-07*t**3

        t1 = (y1-2000.0)/100.0
        t2 = (y2-2000.0)/100.0

        xm1 = xm(t1)
        xn1 = xn(t1)
        xm2 = xm(t2)
        xn2 = xn(t2)

        ra1 = ra/DEGRAD
        dec1 = dec/DEGRAD
        ra2 = ra1+(xm2-xm1)+(xn2-xn1)*sin(ra1)*tan(dec1)
        dec2 = dec1+(xn2-xn1)*cos(ra1)
        ram = (ra1+ra2)/2.
        decm = (dec1+dec2)/2.
        ra2 = ra1+(xm2-xm1)+(xn2-xn1)*sin(ram)*tan(decm)
        dec2 = dec1+(xn2-xn1)*cos(ram)

        ra = ra2*DEGRAD
        dec = dec2*DEGRAD

        return
      end


      subroutine flipi4(iarr,nroe)
!	convert from sun to pc byte order in iarr(nroe)

        integer*4 nroe,iroe
        integer*4 iarr(nroe)
        integer*4 iword
        byte bytes(4)
        byte tempbyte
        equivalence (bytes,iword)

        do iroe = 1,nroe
          iword = iarr(iroe)
          tempbyte = bytes(1)
          bytes(1) = bytes(4)
          bytes(4) = tempbyte
          tempbyte = bytes(2)
          bytes(2) = bytes(3)
          bytes(3) = tempbyte
          iarr(iroe) = iword
        enddo   

        return
      end


      subroutine flipr4(rarr,farr,nroe)
!	farr(nroe) = rarr(nroe) with bytes flipped

        integer*4 nroe,iroe
        real*4 rarr(nroe),farr(nroe)
        real*4 rword
        integer*4 iword
        byte bytes(4)
        byte tempbyte
        equivalence (bytes,iword)
        equivalence (rword,iword)

        do iroe = 1,nroe
          rword = rarr(iroe)
          tempbyte = bytes(1)
          bytes(1) = bytes(4)
          bytes(4) = tempbyte
          tempbyte = bytes(2)
          bytes(2) = bytes(3)
          bytes(3) = tempbyte
          farr(iroe) = rword
        enddo   

        return
      end


      function fnlin(iraw,nsat)
!	convert raw counts to counts/sec with non-linearity correction

        use modes

        integer, parameter :: marr = 16384
        real nodpa,lores,kmirror,krot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        character(60) object,feature,obsmode,instmode
        character(60) objtype,piname,pid,note,weather,warning,version
        real farr(marr)

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
                piname,pid,note,weather,warning,version
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf

        save farr,gainfac,sat,rawmax,zeroval,irawmax

        nfr = nframe*nsum*nwrite
        if(iraw==1) then
!	special case of iraw to initialize parameters and lookup table
          nsat = 0
          gainfac = frtime*gain
          if(index(obsmode,'fowler')>0) gainfac = gainfac*nframe
          if(xnlin==0.) then
            sat = 0.
            irawmax = -satval*nfr
            rawmax = satval
          else
!	make lookup table for fnlin(rawi)
!	rawi, satval, rawmax, farri are in counts/frame
!	if negative, xnlin, ynlin are in counts/sec
!	   and log nlin function is used
!	if positive, they are in counts/frame
!	   and exp nlin function is used
!	jraw is in (counts/4)/frame
            if(xnlin<0.) then
              xnabs = -xnlin*frtime
            else
              xnabs = xnlin
            endif
            if(ynlin<0.) then
              print '(" fnlin is not set up for ynlin<0")'
              stop
            else
              ynabs = ynlin
            endif
            farri = 0.
            draw = 4.
            dfarr = draw
            irawmax = -satval*nfr
            jrawmax = satval/4.
            do jraw = 1,marr
              if(dfarr==0.) then
                farr(jraw) = farrmax
                cycle
              endif
              rawj = jraw*draw
              farri = farri+dfarr
              if(xnlin<0.) then
                rawi = xnabs*log((farri+xnabs)/xnabs)
              elseif(xnlin>0.) then
                rawi = xnabs*(1.-sexp(-farri/xnabs))
              else
                rawi = farri
              endif
              if(znlin>0.) then
                if(ynlin>0.) then
                  gnlin = ynabs*(1.-sexp(-farri/ynabs))
                  rawi = (1.-znlin)*rawi+znlin*gnlin
                else
                  rawi = (1.-znlin)*rawi+znlin*farri
                endif
              endif
              ddfarr = dfarr*(rawj-rawi)/draw
              farri = farri+ddfarr
              dfarr = dfarr+ddfarr
              if((.not.(dfarr>0.)).or.(.not.(rawi>0.))) then
                dfarr = 0.
                farrmax = farri
                if(jraw<jrawmax) then
                  jrawmax = jraw
                  irawmax = -4*jrawmax*nfr
                endif
              endif
              farr(jraw) = farri
            enddo   
            rawmax = farr(jrawmax)
            if(xnlin>0.) then
              sat = xnlin*gain
              irawmax = -0.9*sat*nfr
              rawmax = -sat*alog(0.1)
            else
              sat = 10.*xnabs*gain
              irawmax = -0.9*sat*nfr
              rawmax = -sat*alog(0.1)
            endif
            print '(" non-linearity parameters:",2f10.0,f6.3,i8,2es9.2)', &
              xnlin,ynlin,znlin,-irawmax,rawmax,satval
          endif
!	note: darkval ~ 0 except in rdrst mode
          zeroval = darkval/(frtime*gain)
          if(fowler.and.(darkval>(0.))) zeroval = 0.
          fnlin = 0.
        else
!	convert iraw to corrected counts/frame
          rawi = -iraw/nfr
          if(iraw>1) then
            fnlin = 0.
          elseif(modecard==MNONE) then
            fnlin = rawi
          elseif(iraw<irawmax) then
            fnlin = rawmax/gainfac
            nsat = nsat+1
          elseif(xnlin==0.) then
            fnlin = zeroval+rawi/gainfac
          else
            rawj = rawi/4.
            jraw = rawj
            draw = rawj-jraw
            fnlin = ((1.-draw)*farr(jraw)+draw*farr(jraw+1))/gainfac
          endif
        endif

        return
      end



      subroutine readcard(ierr)
!	read card data, convert to floating point and write to fits

        use ius
        use dims
        use modes

        real cdtr(mc,mx,my),cdmed(mc,mx)
        integer*4 icard(mx),irawmax
        real nodpa,lores,kmirror,krot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical doflip,doflop
        logical badad,badcols,docardfits
        character(32) rawfile,cardfile,lastcard,redfile,sumfile,posfile
        character(48) rawfits,cardfits,redfits,sumfits
        character(4) ired
        character(60) object,feature,obsmode,instmode
        character(60) objtype,piname,pid,note,weather,warning,version

        common /byteflip/ doflip,doflop
        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
                piname,pid,note,weather,warning,version
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /filenames/rawfile,cardfile,lastcard,redfile,sumfile,posfile, &
                rawfits,cardfits,redfits,sumfits,ired
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime

        icbad = 2**30
        ierr = 0

        if(nbitpix/=32) then
          print '(" Bitpix = ",i8," != 32 in ",a24)', nbitpix,cardfile
          write(iupl,'(" Bitpix != 32 in ",a24)') cardfile
          ierr = 1
          return
        endif

        if(.not.((beamtime*abs(gain))>0.)) then
          print '(" Beamtime or gain bad",2e10.3)', beamtime,gain
          write(iupl,'(" Beamtime or gain bad",2e10.3)') beamtime,gain
          ierr = 1
          return
        endif

        if(nc>mc) then
          print '(" Too many card frames for array.",i6, &
            " Read only first",i2)', nc,mc
          nc = mc
        endif

        badad = (date>(7.1001)).and.(date<(7.1019))
        badcols = (frtime<(1.0)) &
          .and.(index(obsmode,'fowler')<=0)

!	open cardfile and cardfits

        open(unit=iurawd,file=cardfile,access='direct',recl=1024, &
            status='old',err=199)

        docardfits = dofits.and.(.not.rdfits) &
          .and.(index(cardfits,'flt.fits')/=1)
        nz = nwrite*nc
        if(docardfits) call rawstart(cardfits,nz,ierr)

!	read card data by rows, write out fits, and convert to floating

        finit = fnlin(1,nsat)
      if(nwrite==1) then

        irec = 0
        lc = 0
        do ic = 1,nc
          do iy = 1,ny
            if(rdhalf.and.((iy<=ny/4).or.(iy>3*ny/4))) then
              do ix = 1,nx
                icard(ix) = -ic
              enddo
            else
              irec = irec+1
              read(unit=iurawd,rec=irec,err=190) icard
              if(doflip) call flipi4(icard,nx)
            endif
            if(docardfits) call fitsidat(iufraw,icard,nx)
            if(badcols) then
              do ix = 1,4
                icard(ix) = -1
              enddo   
            endif
            do ix = 1,nx
              if(abs(icard(ix))>icbad) then
                print '("*** Bad card data.  Is doflip wrong?", &
                  l4,3i4,2i12)', doflip,ic,iy,ix,icard(ix),icbad
                write(iupl,'("!!! Bad card data.  Is doflip wrong?")')
                goto 190
              endif
              cards(ix,iy,ic) = fnlin(icard(ix),nsat)
            enddo   
          enddo   
          if(rdhalf) irec = irec+ny/2
          lc = ic
        enddo   

      elseif(nwrite==2) then

        irec = 0
        lc = 0
        do ic = 1,nc
          do iy = 1,ny
            if(rdhalf.and.((iy<=ny/4).or.(iy>3*ny/4))) then
              do ix = 1,nx
                icard(ix) = -ic
              enddo
            else
              irec = irec+1
              read(unit=iurawd,rec=irec,err=190) icard
              if(doflip) call flipi4(icard,nx)
              if(docardfits) call fitsidat(iufraw,icard,nx)
            endif
            if(badcols) then
              do ix = 1,4
                icard(ix) = -1
              enddo   
            endif
            do ix = 1,nx
              if(iabs(icard(ix))>icbad) then
                print '("*** Bad card data.  Is doflip wrong?", &
                  l4,3i4,2i12)', doflip,ic,iy,ix,icard(ix),icbad
                write(iupl,'("!!! Bad card data.  Is doflip wrong?")')
                goto 190
              endif
              cards(ix,iy,ic) = fnlin(icard(ix),nsat)
            enddo   
          enddo   
          if(rdhalf) irec = irec+ny/2

          do iy = 1,ny
            if(rdhalf.and.((iy<=ny/4).or.(iy>3*ny/4))) then
              do ix = 1,nx
                icard(ix) = -ic
              enddo
            else
              irec = irec+1
              read(unit=iurawd,rec=irec,err=190) icard
              if(doflip) call flipi4(icard,nx)
            endif
            if(docardfits) call fitsidat(iufraw,icard,256)
            if(badcols) then
              do ix = 1,4
                icard(ix) = -1
              enddo   
            endif
            do ix = 1,nx
              cardi = fnlin(icard(ix),nsat)
!	try taking the fainter to avoid spikes
              cards(ix,iy,ic) = 2.*min(cards(ix,iy,ic),cardi)
            enddo   
          enddo   
          if(rdhalf) irec = irec+ny/2

          lc = ic
        enddo   

      else

!	average multiple cards.  median might be better
!	might also check for saturated pixels

        do ic = 1,nc
        do iy = 1,ny
        do ix = 1,nx
          cards(ix,iy,ic) = 0.
        enddo   
        enddo   
        enddo   

        irec = 0
        lc = 0
        do ic = 1,nc
          do iw = 1,nwrite
            ifr = (ic-1)*nwrite+iw
            do iy = 1,ny
              irec = irec+1
              read(unit=iurawd,rec=irec,err=190) icard
              if(doflip) call flipi4(icard,nx)
              if(docardfits) call fitsidat(iufraw,icard,256)
              if(badcols) then
                do ix = 1,4
                  icard(ix) = -1
                enddo   
              endif
              do ix = 1,nx
                if(iabs(icard(ix))>icbad) then
                  print '("*** Bad card data.  Is doflip wrong?", &
                    l4,i4,2i12)', doflip,ix,icard(ix),icbad
                  write(iupl,'("!!! Bad card data.  Is doflip wrong?")')
                  goto 190
                endif
                cards(ix,iy,ic) = cards(ix,iy,ic)+fnlin(icard(ix),nsat)
              enddo   
            enddo   
          enddo   
          lc = ic
        enddo   

      endif
      if(docardfits) then
        call fitsipad(iufraw)
        close(unit=iufraw)
      elseif(rdfits) then
        close(unit=iufraw)
      endif

        if(nsat>0) &
          print '(i8," pixels exceed non-linearity max")', nsat

        if(fowler.and.(darkval>=0.)) then
          xlc = lc
          if((darkval==0.).or.(darkval>xlc)) then
            icmn = 0
            summn = 1.e+32
            do ic = 1,lc
              sumc = 0.
              do iy = 1,ny
                do ix = 1,nx
                  sumc = sumc+cards(ix,iy,ic)
                enddo   
              enddo   
              if(sumc<summn) then
                icmn = ic
                summn = sumc
              endif
            enddo   
          else
            icmn = nint(darkval)
          endif
          if(icmn>0) then
            do ix = 1,nx
              sumc = 0.
              do iy = 1,4
                sumc = sumc+cards(ix,iy,icmn)
              enddo   
              sumc = sumc/4.
              do iy = 1,ny
                do ic = 1,lc
                  cards(ix,iy,ic) = cards(ix,iy,ic)-sumc
                enddo   
              enddo   
            enddo   
          endif
        endif

        ierr = 0
        if(docardfits) call fitsclos(iufraw)
        close(unit=iurawd)
        return

  190   if((lc>=2).or.((lc==1).and.(modecard==MBLKOBJ))) then
          print '(" Read only ",i1,"/",i1," flat frames ",i6)',lc,nc,irec
          if(lc<nc) ierr = 1
          nc = lc
        else
          print '(" Error reading data from ",a16)', cardfile
          write(iupl,'(" Error reading data from ",a16)') cardfile
          ierr = 9
        endif
        if(docardfits) call fitsclos(iufraw)
        close(unit=iurawd)
        return

  199   print '(" Error opening file ",a16)', cardfile
        write(iupl,'(" Error opening file ",a16)') cardfile
        ierr = 99
        close(unit=iurawd)
        return

      end


      subroutine marscard(ab,bb,ierr)
!	read marscard data and convert to normal card format

        use ius
        use dims
        use modes

        real ab(mx,my,mp),bb(mx,my,mp)
        integer*4 irawmax
        real nodpa,lores,kmirror,krot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical dokeep
        character(32) rawfile,cardfile,lastcard,redfile,sumfile,posfile
        character(48) rawfits,cardfits,redfits,sumfits
        character(4) ired

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /filenames/rawfile,cardfile,lastcard,redfile,sumfile,posfile, &
                rawfits,cardfits,redfits,sumfits,ired
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime
        common /iwin/ iwin1,iwin2

        if(nbitpix/=32) then
          print '(" Bitpix = ",i8," != 32 in ",a24)', nbitpix,cardfile
          write(iupl,'(" Bitpix != 32 in ",a24)') cardfile
          ierr = 1
          return
        endif

        if(.not.((beamtime*abs(gain))>0.)) then
          print '(" Beamtime or gain bad",2e10.3)', beamtime,gain
          write(iupl,'(" Beamtime or gain bad",2e10.3)') beamtime,gain
          ierr = 1
          return
        endif

        do ic = 1,mc
          do iy = 1,ny
            do ix = 1,nx
              cards(ix,iy,ic) = 0.
            enddo
          enddo
        enddo

        dokeep = dofits
        dofits = .false.

        ierr = 0
        jerr = -1
      if((date>(3.05)).and.(date<(4.00))) then
!	marsflat contains repeated black, shiny sequence
!	use readraw to read black and shiny into bb and ab
!	then put median filtered frames into cards(,,1 and ,,3)

        nodon = .true.
        np = -4
        call readraw(cardfile,ab,bb,np,ierr)
        iq = 2
        call median(cards(1,1,1),bb,mx,my,np,iq,ierr)
        call median(cards(1,1,3),ab,mx,my,np,iq,ierr)
          print '(" Plotting medians")'
          call grey(mx,nx,ny,cards(1,1,1),0.,0.,iwin1)
          call waitasec(pause,.true.,cards(1,1,1),mx,nx,ny)
          call grey(mx,nx,ny,cards(1,1,3),0.,0.,iwin1)
          call waitasec(pause,.true.,cards(1,1,3),mx,nx,ny)
        if(verbose) then
          print '(" Plotting marscard sequence")'
          do ip = 1,np
            call grey(mx,nx,ny,bb(1,1,ip),0.,0.,iwin1)
            call waitasec(pause,.true.,bb(1,1,ip),mx,nx,ny)
            call grey(mx,nx,ny,ab(1,1,ip),0.,0.,iwin1)
            call waitasec(pause,.true.,ab(1,1,ip),mx,nx,ny)
          enddo
          print '(" Plotting medians")'
          call grey(mx,nx,ny,cards(1,1,1),0.,0.,iwin1)
          call waitasec(pause,.true.,cards(1,1,1),mx,nx,ny)
          call grey(mx,nx,ny,cards(1,1,3),0.,0.,iwin1)
          call waitasec(pause,.true.,cards(1,1,3),mx,nx,ny)
        endif

!	use readraw to read sky from rawfile into ab
!	then put quartile filtered frame into cards(,,2)

        call readhdr(.true.,jerr)
        np = -3
        call readraw(rawfile,ab,bb,np,ierr)
        iq = 1
        call median(cards(1,1,2),ab,mx,my,np,iq,ierr)
        if(verbose) then
          print '(" Plotting sky median")'
          call grey(mx,nx,ny,cards(1,1,2),0.,0.,iwin1)
          call waitasec(pause,.true.,cards(1,1,2),mx,nx,ny)
        endif
        nc = 3

      elseif((date>(4.00)).and.(date<(5.03))) then
!	manual flats containing just black
!	use readraw to read black into ab
!	then put median filtered frame into cards(,,1)

        nodon = .true.
        np = -3
        call readraw(cardfile,ab,bb,np,ierr)
        iq = 2
        call median(cards(1,1,1),ab,mx,my,np,iq,ierr)
        do iy = 1,ny
          do ix = 1,nx
            cards(ix,iy,3) = cards(ix,iy,1)
          enddo
        enddo

!	use readraw to read sky from rawfile into ab
!	then put quartile filtered frame into cards(,,2)

        call readhdr(.true.,jerr)
        np = -3
        call readraw(rawfile,ab,bb,np,ierr)
        iq = 1
        call median(cards(1,1,2),ab,mx,my,np,iq,ierr)
        nc = 2

      elseif((date>(5.07)).and.(date<(6.00))) then
!	8 (or np/2) blacks followed by 8 skys
!	use readraw to read into ab
!	then put median filtered frames 1-8 into cards(,,1)
!	and put median filtered frames 9-16 into cards(,,2)

        nodon = .true.
        np = -3
        call readraw(cardfile,ab,bb,np,ierr)
        nc = np/2
        iq = 2
        call median(cards(1,1,1),ab(1,1,1),mx,my,nc,iq,ierr)
        call median(cards(1,1,2),ab(1,1,nc+1),mx,my,nc,iq,ierr)
        do iy = 1,ny
        do ix = 1,nx
          cards(ix,iy,3) = cards(ix,iy,1)
        enddo
        enddo
        nc = 2

      else
        print '("*** I don''t know how to process marsflats for ", &
          f8.4)', date
        ierr = 9
      endif

        dofits = dokeep

        return
      end


      function xmed(arr,nx,iq,ierr)
!	return median or quartile of 1-D array

        real arr(nx),arrx(16384)

        nxq = nx
        arrx(1) = arr(1)
        if(arrx(1)==0.) nxq = nxq-1
        do ix = 2,nx
          arri = arr(ix)
          if(arri==0.) then
            nxq = nxq-1
            arrx(ix) = arri
            cycle
          endif
          do jx = ix,2,-1
            arrj = arrx(jx-1)
            if((arrj==0.).or.(arri<arrj)) then
              arrx(jx) = arrj
            else
              arrx(jx) = arri
              goto 100
            endif
          enddo
          arrx(1) = arri
  100   enddo
        if(iq==0) then
          xmed = arrx(1)
        elseif((iq<0).and.(-iq<=nxq)) then
          xmed = arrx(-iq)
        elseif((iq>4).and.(iq-5<nzq)) then
          xmed = arrx(nzq-(iq-5))
        elseif(iq==4) then
          ixq = nxq-(nxq/64)
          xmed = arrx(ixq)
        else
          ixq = (iq*nxq+3)/4
          if((mod(nxq,4)==0).or.((iq==2).and.(mod(nxq,2)==0))) then
            xmed = (arrx(ixq)+arrx(ixq+1))/2.
          else
            xmed = arrx(ixq)
          endif
        endif

      end


      subroutine median(arr2,arr3,nx,ny,nz,iq,ierr)
!	collapse arr3 in z by median or quartile filtering

        use dims

        real arr2(nx,ny),arr3(nx,ny,nz)
        real arrz(mp)

        do iy = 1,ny
          do ix = 1,nx
            nzq = nz
            arrz(1) = arr3(ix,iy,1)
            if(arrz(1)==0.) nzq = nzq-1
            do iz = 2,nz
              arri = arr3(ix,iy,iz)
              if(arri==0.) then
                nzq = nzq-1
                arrz(iz) = arri
                goto 100
              endif
              do jz = iz,2,-1
                arrj = arrz(jz-1)
                if((arrj==0.).or.(arri<arrj)) then
                  arrz(jz) = arrj
                else
                  arrz(jz) = arri
                  goto 100
                endif
              enddo   
              arrz(1) = arri
  100       enddo   
            if(iq==0) then
              arr2(ix,iy) = arrz(1)
            elseif((iq<0).and.(-iq<=nzq)) then
              arr2(ix,iy) = arrz(-iq)
            elseif(iq==4) then
              izq = nzq-(nzq/64)
              arr2(ix,iy) = arrz(izq)
            elseif((iq>4).and.(iq-5<nzq)) then
              arr2(ix,iy) = arrz(nzq-(iq-5))
            else
              izq = (iq*nzq+3)/4
              if((mod(nzq,4)==0).or.(iq==2).and.(mod(nzq,2)==0)) then
                arr2(ix,iy) = (arrz(izq)+arrz(izq+1))/2.
              else
                arr2(ix,iy) = arrz(izq)
              endif
            endif
          enddo   
        enddo   

        return
      end



      subroutine objcard(rawfile,ab,bb,np,ierr)
!	read data with flats inserted in object file if(ierr==4)
!	or just get sky from object file and insert into flat array

        use ius
        use dims
        use modes

        character(32) rawfile
        real ab(mx,my,mp),bb(mx,my,mp)
        integer*4 irawmax
        real nodpa,lores,kmirror,krot
        logical baddata,doplot
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime
        common /iwin/ iwin1,iwin2

        doplot = (iwin1>0)
        iopt = ierr

        if(nbitpix/=32) then
          print '(" Bitpix = ",i8," != 32 in ",a24)', nbitpix,cardfile
          write(iupl,'(" Bitpix != 32 in ",a24)') cardfile
          ierr = 1
          return
        endif

        if(.not.((beamtime*abs(gain))>0.)) then
          print '(" Beamtime or gain bad",2e10.3)', beamtime,gain
          write(iupl,'(" Beamtime or gain bad",2e10.3)') beamtime,gain
          ierr = 1
          return
        endif

!	used readraw to read sky from rawfile into ab
!	if iopt = 4 put brightest frame into cards(,,1) and cards(,,3)
!	if iopt = 2-4 put first quadrant frame into cards(,,2) and cards(,,4)

        if((modeobs==MSCAN).and.(intsky(1)>0) &
          .and.(intsky(2)>=intsky(1))) then
          ip1 = intsky(1)
          ip2 = intsky(2)
          np = ip2-ip1+1
        else
          if((modeobs==MNOD).or.(modeobs==MCHOP)) np = 2*np
          ip1 = 1
          ip2 = np
        endif
        print '(" Read",i4," sky-flat frames from ",a16)', np,rawfile
        if(iopt==4) then
          iq = 5
          call median(cards(1,1,1),ab,mx,my,np,iq,ierr)
          do iy = 1,ny
            do ix = 1,nx
              cards(ix,iy,3) = cards(ix,iy,1)
            enddo
          enddo
        endif
        if(iopt==1) then
          iq = 2
          call median(cards(1,1,2),ab,mx,my,np,iq,ierr)
          if(verbose) then
            print '(" Plotting medians(sky)")'
            call grey(mx,nx,ny,cards(1,1,2),0.,0.,iwin1)
            call waitasec(pause,.true.,cards(1,1,2),mx,nx,ny)
          endif
        else
          iq = 1
          call median(cards(1,1,2),ab,mx,my,np,iq,ierr)
          if(verbose) then
            print '(" Plotting medians(sky)")'
            call grey(mx,nx,ny,cards(1,1,2),0.,0.,iwin1)
            call waitasec(pause,.true.,cards(1,1,2),mx,nx,ny)
          endif
        endif
        do iy = 1,ny
          do ix = 1,nx
            cards(ix,iy,4) = cards(ix,iy,2)
          enddo
        enddo
        if(iopt==3) then
!	use old blk with sky from object file
          do iy = 1,ny
            do ix = 1,nx
              cards(ix,iy,1) = cards(ix,iy,3)
            enddo
          enddo
        elseif(iopt==2) then
!	flat field with black-objsky
          do iy = 1,ny
            do ix = 1,nx
              cards(ix,iy,3) = cards(ix,iy,1)
            enddo
          enddo
        endif
        nc = 2

        return
      end


      subroutine readraw(rdfile,ab,bb,np,ierr)
!	try to read nnod raw data pairs, and return np = number of pairs read
!	if(dofits.and.(.not.rdfits)) copy rdfile to rawfits

        use ius
        use dims
        use modes

        character(32) rdfile
        real ab(mx,my,mp),bb(mx,my,mp),vali(mx)
        character(60) object,feature,obsmode,instmode
        character(60) objtype,piname,pid,note,weather,warning,version
        real nodpa,lores,kmirror,krot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        character(32) rawfile,cardfile,lastcard,redfile,sumfile,posfile
        character(48) rawfits,cardfits,redfits,sumfits
        character(4) ired
        character(80) line
        integer*4 iraw(mx),irawmax
        logical doflip,doflop
        logical alla,lastscan,baab
        logical badad,badcols
        logical dorawfits,rdrawfits

        common /byteflip/ doflip,doflop
        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
                piname,pid,note,weather,warning,version
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /filenames/rawfile,cardfile,lastcard,redfile,sumfile,posfile, &
                rawfits,cardfits,redfits,sumfits,ired

        save irec,nzp
        rdrawfits = rdfits.or.(index(rdfile,'.fits')>1)
        dorawfits = (index(rawfits,'raw.fits')/=1).and.(.not.rdrawfits)

        if(np<0) then
          irec = 0
          nzp = 0
        endif
        alla = ((modeobs==MSTARE).or.(modeobs==MFLAT) &
             .or.(modeobs==MSCAN).or.(np==(-3)))
        baab = abba.and.(.not.alla).and.(nwrite==2)
        if(np==-2) then
!	reading object file to get sky for flat
!          alla = .false.
          alla = .true.
          if((modobs==MNOD).or.(modeobs==MMAP)) nnod = 2*nnod
        endif
        if((np==-3).and.((modeobs==MNOD).or.(modeobs==MMAP))) &
          nnod = 2*nnod
        if(np==-4) alla = .false.
        np = 0
        if((modeobs==MCHOP).or.(modeobs==MNOD).or.(modeobs==MCHOPNOD)) then
          nz = 2*nnod
        elseif(modeobs==MSCAN) then
          nz = nnod*nscan
        elseif(modeobs==MMAP) then
          nz = 2*nnod*nscan
        else
          nz = nnod
        endif

        if(nbitpix/=32) then
          print '(" Bitpix = ",i8," != 32 in ",a24)', nbitpix,rdfile
          write(iupl,'(" Bitpix != 32 in ",a24)') rdfile
          ierr = 1
          return
        endif

        if(.not.((beamtime*abs(gain))>0.)) then
          print '(" Beamtime or gain bad",2e10.3)', beamtime,gain
          write(iupl,'(" Beamtime or gain bad",2e10.3)') beamtime,gain
          ierr = 1
          return
        endif

        badad = (date>(7.1001)).and.(date<(7.1019))
        badcols = (frtime<(1.0)) &
          .and.(index(obsmode,'fowler')<=0)

        if(modeobs==MSCAN) then
          lastscan = (irec>=(nnod-nsky)*nwrite*ny*(nscan-1))
        elseif(modeobs==MMAP) then
          lastscan = (irec>=2*nnod*nwrite*ny*(nscan-1))
        else
          lastscan = .true.
        endif

        if(ierr==(-1)) then
!	skip this frame or scan
          if(modeobs==MSCAN) then
            irec = irec+(nnod-nsky)*nwrite*ny
          elseif(modeobs==MMAP) then
            irec = irec+2*nnod*nwrite*ny
          elseif((modeobs==MSTARE).or.(modeobs==MFLAT)) then
            print '(" How did I get here?")'
            irec = irec+1
          else
            print '(" How did I get here?")'
            irec = irec+2
          endif
          ierr = 0
          goto 900
        endif

!	open rdfile and start rawfits file

        lrec = 4*nx
        open(unit=iurawd,file=rdfile, &
          access='direct',recl=lrec,status='old',err=199)
        if(dorawfits.and.(nzp==0)) then
          call rawstart(rawfits,nz,ierr)
        endif
        if(ierr/=0) then
          print '("*** Error starting ",a36," Won''t store data")', rawfits
          dorawfits = .false.
          ierr = 0
        endif

!	read raw data by rows, convert to floating
!	and (try to) correct for detector non-linearity

        finit = fnlin(1,nsat)

      if(nwrite==1) then

        if(modeobs==MSCAN) then
          ifr = nnod*(irec/ny)/(nnod-nsky)
        else
          ifr = irec/ny
        endif
        do in = 1,nnod
          ifr = ifr+1
          do iy = 1,ny
            if(rdhalf.and.((iy<=ny/4).or.(iy>3*ny/4))) then
              do ix = 1,nx
                iraw(ix) = -in
              enddo
            else
              irec = irec+1
              read(unit=iurawd,rec=irec,err=190) iraw
              if(doflip) call flipi4(iraw,nx)
            endif
            if(dorawfits) call fitsidat(iufraw,iraw,nx)
            if(badcols) then
              do ix = 1,4
                iraw(ix) = -1
              enddo   
            endif
            if(badad) then
              do ix = 3,nx,4
                iraw(ix) = (iraw(ix-1)+iraw(ix+1))/2
              enddo   
            endif
            if(alla.or.(date<(1.00))) then
              do ix = 1,nx
                ab(ix,iy,in) = fnlin(iraw(ix),nsat)
              enddo   
            else
!	first frame is in B beam in 2001 and after
              do ix = 1,nx
                bb(ix,iy,in) = fnlin(iraw(ix),nsat)
              enddo   
            endif
          enddo   
          if(rdhalf) irec = irec+ny/2

          if((modeobs==MSCAN).and.(dosubsky.or.(intsky(1)>=0))) then
            do iy = 1,ny
              do ix = 1,nx
                bb(ix,iy,in) = ab(ix,iy,1)+1.0/gain
              enddo   
            enddo   
          elseif(alla) then
            do iy = 1,ny
              do ix = 1,nx
                bb(ix,iy,in) = 0.
              enddo   
            enddo   
          else
            ifr = ifr+1
            do iy = 1,ny
              if(rdhalf.and.((iy<=ny/4).or.(iy>3*ny/4))) then
                do ix = 1,nx
                  iraw(ix) = -in
                enddo
              else
                irec = irec+1
                read(unit=iurawd,rec=irec,err=190) iraw
                if(doflip) call flipi4(iraw,nx)
              endif
              if(dorawfits) call fitsidat(iufraw,iraw,nx)
              if(badcols) then
                do ix = 1,4
                  iraw(ix) = -1
                enddo   
              endif
              if(badad) then
                do ix = 3,nx,4
                  iraw(ix) = (iraw(ix-1)+iraw(ix+1))/2
                enddo   
              endif
              if(date<(1.00)) then
                do ix = 1,nx
                  bb(ix,iy,in) = fnlin(iraw(ix),nsat)
                enddo   
              else
                do ix = 1,nx
                  ab(ix,iy,in) = fnlin(iraw(ix),nsat)
                enddo   
              endif
            enddo   
            if(rdhalf) irec = irec+ny/2
          endif

          np = np+1
        enddo   

      elseif(baab) then
        print '(" Treating nwrite=2 data as alternating ba, ab pairs")'

        nwrite = nwrite/2
        nnod = 2*nnod-1
        print '(" Changing nwrite, nnod to ",2i4)', nwrite,nnod

!	skip first (B) frame
        irec = irec+ny
        do in = 1,nnod
          ifr = ifr+1
          do iy = 1,ny
            irec = irec+1
            read(unit=iurawd,rec=irec,err=190) iraw
            if(doflip) call flipi4(iraw,nx)
            if(dorawfits) call fitsidat(iufraw,iraw,nx)
            if(badcols) then
              do ix = 1,4
                iraw(ix) = -1
              enddo   
            endif
            if(mod(in,2)==1) then
              do ix = 1,nx
                bb(ix,iy,in) = fnlin(iraw(ix),nsat)
              enddo   
            else
              do ix = 1,nx
                ab(ix,iy,in) = fnlin(iraw(ix),nsat)
              enddo   
            endif
          enddo   

          ifr = ifr+1
          do iy = 1,ny
            irec = irec+1
            read(unit=iurawd,rec=irec,err=190) iraw
            if(doflip) call flipi4(iraw,nx)
            if(dorawfits) call fitsidat(iufraw,iraw,nx)
            if(badcols) then
              do ix = 1,4
                iraw(ix) = -1
              enddo   
            endif
            if(mod(in,2)==0) then
              do ix = 1,nx
                bb(ix,iy,in) = fnlin(iraw(ix),nsat)
              enddo   
            else
              do ix = 1,nx
                ab(ix,iy,in) = fnlin(iraw(ix),nsat)
              enddo   
            endif
          enddo   

          np = np+1
        enddo   

      else
!	nwrite > 1 and not treating data as abba sequences

        do in = 1,nnod
          do iy = 1,ny
            do ix = 1,nx
              ab(ix,iy,in) = 0.
              bb(ix,iy,in) = 0.
            enddo   
          enddo   

          if(modeobs==MSCAN) then
            ifr = nnod*(irec/ny)/(nnod-nsky)
          else
            ifr = irec/ny
          endif
          do iw = 1,nwrite
            ifr = ifr+1
            do iy = 1,ny
              irec = irec+1
              read(unit=iurawd,rec=irec,err=190) iraw
              if(doflip) call flipi4(iraw,nx)
              if(dorawfits) call fitsidat(iufraw,iraw,nx)
              if(badcols) then
                do ix = 1,4
                  iraw(ix) = -1
                enddo   
              endif
              if(alla.or.(date<(1.00))) then
                do ix = 1,nx
                  ab(ix,iy,in) = ab(ix,iy,in)+fnlin(iraw(ix),nsat)
                enddo   
              else
                do ix = 1,nx
                  bb(ix,iy,in) = bb(ix,iy,in)+fnlin(iraw(ix),nsat)
                enddo   
              endif
            enddo   
          enddo   

          if((modeobs==MSCAN).and.(dosubsky.or.(intsky(1)>=0))) then
            do iy = 1,ny
              do ix = 1,nx
                bb(ix,iy,in) = ab(ix,iy,1)+1.0/gain
              enddo   
            enddo   
          elseif(.not.alla) then
            ifr = ifr+1
            do iw = 1,nwrite
              ifr = ifr+1
              do iy = 1,ny
                irec = irec+1
                read(unit=iurawd,rec=irec,err=190) iraw
                if(doflip) call flipi4(iraw,nx)
                if(dorawfits) call fitsidat(iufraw,iraw,nx)
                if(badcols) then
                  do ix = 1,4
                    iraw(ix) = -1
                  enddo   
                endif
                if(date<(1.00)) then
                  do ix = 1,nx
                    bb(ix,iy,in) = bb(ix,iy,in)+fnlin(iraw(ix),nsat)
                  enddo   
                else
                  do ix = 1,nx
                    ab(ix,iy,in) = ab(ix,iy,in)+fnlin(iraw(ix),nsat)
                  enddo   
                endif
              enddo   
            enddo   
          endif

          np = np+1
        enddo   

      endif

  400   if(.not.(alla.or.nodon)) then
          if(date<(1.00)) then
            do in = np,2,-1
              in2 = in-1
              do iy = 1,ny
                do ix = 1,nx
                  bb(ix,iy,in) = &
                    (bb(ix,iy,in)+bb(ix,iy,in2))/2.
                enddo   
              enddo   
            enddo   
          elseif(modeobs==MMAP) then
            do in = 1,np-1
              in2 = in+1
              do iy = 1,ny
                do ix = 1,nx
                  bb(ix,iy,in) = &
                    (bb(ix,iy,in)+bb(ix,iy,in2))/2.
                enddo   
              enddo   
            enddo   
          endif
        endif

        print '(" Read ",i3," pairs from ",a16)', np,rdfile
        write(iupl,'(" Read ",i3," pairs from ",a16)') np,rdfile
        if(nsat>0) &
          print '(i8," pixels exceed non-linearity max")', nsat
        if(dorawfits) then
          if((modeobs==MCHOP).or.(modeobs==MNOD) &
            .or.(modeobs==MCHOPNOD).or.(modeobs==MMAP)) then
            nzp = nzp+2*np
          else
            nzp = nzp+np
          endif
          if(lastscan.and.(nzp<nz)) then
            write(line,'("NAXIS3  = ",i20, &
              "  / edited due to incomplete data file")') nzp
            call fithedit(iufraw,line)
          endif
        endif

        ierr = 0
        if((np<1).or.((np<nnod).and. &
          ((modeobs==MSCAN).or.(modeobs==MMAP)))) then
          lastscan = .true.
          ierr = 8
        elseif((np==nnod).and. &
          ((modeobs==MSCAN).or.(modeobs==MMAP))) then
          irec1 = irec+1
          read(unit=iurawd,rec=irec1,err=900) iraw
          print '(" Did not read to end of file ",a16)', rdfile
          write(iupl,'(" Did not read to end of file ",a16)') rdfile
        endif
        if(modeobs==MSCAN) irec = irec-nwrite*ny*nsky
        goto 900

  190   if((modeobs/=MSCAN).and.(modeobs/=MMAP)) then
          lastscan = .true.
          goto 400
        else
          print '("*** Error reading from ",a16)', rdfile
          write(iupl,'("!!! Error reading from ",a16)') rdfile
          write(iuredh,'("warning = Incomplete data file")')
          ierr = 9
          if(dorawfits.and.(nzp<nz)) then
            write(line,'("NAXIS3  = ",i20, &
              "  / edited due to incomplete data file")') nzp
            call fithedit(iufraw,line)
          endif
          lastscan = .true.
!          close(unit=iuredh)
          goto 900
        endif

  199   print '(" Error opening ",a16)', rdfile
        write(iupl,'("!!! Error opening ",a16)') rdfile
        ierr = 9

  900   if(dorawfits.and.(((modeobs/=MSCAN).and.(modeobs/=MMAP)) &
          .or.lastscan)) then
          call fitsipad(iufraw)
          call fitsclos(iufraw)
        endif
        close(unit=iurawd)
        return

      end


      subroutine makeflat(flat,ierr)
!	use card sequence to make a calibrated flat
!	also run testtort on black or black-shiny
!	to determine distortion parameters and illum

        use ius
        use dims
        use modes

        real flat(mx,my)
        real diff(mx,my),std(mx,my)
        real csum(mc),cmin(mc),cmed(mc)
        real nodpa,lores,kmirror,krot,lskrot
        character(60) object,feature,obsmode,instmode
        character(60) objtype,piname,pid,note,weather,warning,version
        character(8) yn
        logical good
        logical baddata,doplot
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
                piname,pid,note,weather,warning,version
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
        common /tort0/ hrfl0,xdfl0,xdmrdgr,xdlrdgr,xdmr0,xdlr0,fred0, &
                xdkrot,lskrot
        common /iwin/ iwin1,iwin2
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime

        doplot = (iwin1>0)

        if((modecard==MNONE).and.((modetort==MNONE) &
          .or.(modetort==MOBJ).or.(modetort==MCELL) &
          .or.(modetort==MOLD))) then
          nc = 4
          do iy = 1,ny
            do ix = 1,nx
              do ic = 1,nc
                cards(ix,iy,ic) = 1.
              enddo   
              flat(ix,iy) = 1.
              illum(ix,iy) = 1
            enddo   
          enddo   
          return
        endif

        if(modeinst==MCAMERA) then
          modeold = modecard
          if(modecard/=MSKY) then
            print '(" Setting flatmode = shiny for camera mode")'
            modecard = MSHINY
          endif
        endif

        if(nc<3) then
          if(modecard==MSHINY) then
            print '(" Warning: shiny not read.", &
              "  Changing cardmode to blk.")'
            modecard = MBLK
          elseif(modecard==BLKSHINY) then
            print '(" Warning: shiny not read.", &
              "  Changing cardmode to blksky.")'
            modecard = MBLKSKY
          endif
        endif

        if((modecard==MBLKSKY).or.(modecard==MBLKOBJ)) then
          if(dtsky==0.) then
            skyscl = 1.0
          else
            skyscl = bnu(waveno0,temp,ierr)/bnu(waveno0,temp-dtsky,ierr)
          endif
        else
          skyscl = 1.0
        endif

        if(verbose) then
         do ic = 1,nc
          print '(" Plotting raw cards(",i1,")")', ic
          if(pause.and.ask) then
            call grey(mx,nx,ny,cards(1,1,ic),-1.,-1.,iwin1)
          else
            call grey(mx,nx,ny,cards(1,1,ic),0.,0.,iwin1)
          endif
          call waitasec(pause,.true.,cards(1,1,ic),mx,nx,ny)
         enddo   
        endif

!	calculate means to find black (brightest) frame

        do ic = 1,nc
          sum = 0.
          cimin = darkcard
          do iy = 1,ny
            do ix = 1,nx
              sum = sum+cards(ix,iy,ic)
              if(cards(ix,iy,ic)<cimin) cimin = cards(ix,iy,ic)
            enddo   
          enddo
          if((sum/=sum).or.(sum==0.)) then
            print '("*** Bad data in flat frame ",i1,f9.2)', ic,sum
            write(iupl,'("!!! Bad data in flat frame ",i1,f9.2)') ic,sum
            if(ic==3) then
              if((modecard==MSHINY).or.(modecard==MBLKSHINY)) then
                print '("  Changing cardmode to blksky.")'
                modecard = MBLKSKY
              else
                print '(" But I don''t care")'
              endif
            else
              print '(" Setting flat = 1")'
              call setone(cards,nc,ierr)
              call setone(flat,1,ierr)
              ierr = 1
              return
            endif
          endif
          csum(ic) = sum
          cmin(ic) = cimin
        enddo   

!	use quartile filtered values instead?
        do ic = 1,nc
          cmed(ic) = quartf(cards(1,1,ic),mx*my,3,ierr)
        enddo

        if((modeblk>0).and.(modeblk<=nc)) then
          icblk = modeblk
        else
          cmax = csum(1)
          icblk = 1
          do ic = 2,nc
            if(csum(ic)>cmax) then
              cmax = csum(ic)
              icblk = ic
            endif
          enddo   
        endif
!        cblk = csum(icblk)/(nx*ny)
        cblk = cmed(icblk)
        eperadu = 100.
        if(temp==273.16) print '(" Warning: temp not read.  Using 0 C")'
        pnut = pnu(waveno0,temp,ierr)
        aomega = pixelwd**2*(3.14159/4.)/36.
        if((modeinst==MHIMED).or.(modeinst==MHILOW)) then
          dwno = waveno0*slitval(slit,date)/(2.*abs(hrr)*hrfl0)
          rqe = (cblk*eperadu)/(pnut*aomega*dwno)
          print '(" Median RQE for black = ",es9.2, &
            " (assuming",f6.0," e/ADU)")', rqe,eperadu
        elseif((modeinst==MMED).or.(modeinst==MLOW)) then
          dwno = waveno0*slitval(slit,date)/(2.*xdr*xdfl0)
          rqe = (cblk*eperadu)/(pnut*aomega*dwno)
          print '(" Median RQE for black = ",es9.2, &
            " (assuming",f6.0," e/ADU)")', rqe,eperadu
        endif

        if(icblk==1) then
          icsky = 2
          icshiny = 3
          if((nc==4).and.(csum(4)<csum(2))) icsky = 4
        elseif(icblk==2) then
          icsky = 1
          if((nc==4).and.(csum(3)<csum(1))) icsky = 3
          icshiny = 4
        else
          icsky = icblk-1
          icshiny = icblk-2
        endif
        icsky2 = icsky+2
        if(icsky>2) icsky2 = icsky-2
        if(icsky2>nc) icsky2 = icsky
        if((icblk/=1).and.(icblk/=modeblk)) then
          print '(" ** Flat frame #",i1," is brightest")', icblk
          write(iupl,'(" ** Flat frame #",i1," is brightest")') icblk
        endif
        if(icshiny>nc) icshiny = icblk
        cblk = cblk*frtime*gain
        csky = csum(icsky)*frtime*gain/(nx*ny)
        print '(" Mean counts on black, sky =",2f7.0)', cblk,csky
        if(doplot.and.(cblk<satval/4.).and.(csky<1000.)) &
          print '("*** May want to increase frame time")'

        if(cmax<0.) then
          print '("*** All cards < 0.  Try subtracting min")'
          write(iupl,'("!!! All cards < 0.  Try subtracting min")')
          if(modecard==MSKY) then
            czero = cmin(icsky)
          elseif(modecard==MSHINY) then
            czero = cmin(icshiny)
          else
            czero = cmin(icblk)
          endif
        else
          czero = 0.
        endif

!	cmax = satval/(frtime*gain)
!	if((xnlin>0.).and.(xnlin<=1.)) cmax = cmax*16.**(1./xnlin)
        irawmax = -nframe*nsum*nwrite*int(satval)
        cmax = 0.99*fnlin(irawmax,msat)
        eperadu = 100.*beamtime
        readvar = (30./eperadu)**2
        nsat = 0
        sumbl = 0.
        sumsh = 0.
        sumsk = 0.
        sumsl = 0.
        sumbl2 = 0.
        sumsh2 = 0.
        sumsk2 = 0.
        sumsl2 = 0.
        sumblsh = 0.
        sumsksl = 0.
        do ix = 1,nx
          do iy = 1,ny
            sumbl = sumbl+cards(ix,iy,icblk)
            sumsh = sumsh+cards(ix,iy,icshiny)
            sumsk = sumsk+cards(ix,iy,icsky)
            sumsl = sumsl+cards(ix,iy,icsky2)
            sumbl2 = sumbl2+cards(ix,iy,icblk)**2
            sumsh2 = sumsh2+cards(ix,iy,icshiny)**2
            sumsk2 = sumsk2+cards(ix,iy,icsky)**2
            sumsl2 = sumsl2+cards(ix,iy,icsky2)**2
            sumblsh = sumblsh+cards(ix,iy,icblk)*cards(ix,iy,icshiny)
            sumsksl = sumsksl+cards(ix,iy,icsky)*cards(ix,iy,icsky2)
            if((modecard==MBLK).or.(modecard==MNONE)) then
              if((cards(ix,iy,icblk)>cmax).and.(cmax>0.)) then
                good(ix,iy) = .false.
                nsat = nsat+1
              endif
              cardd = cards(ix,iy,icblk)
              card1 = cardd
              if(cardd>0.) then
                card2 = 1.-skyscl*cards(ix,iy,icsky)/cardd
              else
                card2 = 0.
              endif
              cardstd = sqrt(abs(2.*cardd/eperadu)+readvar)
            elseif((modecard==MSHINY).and.(icshiny/=icblk)) then
              if((cards(ix,iy,icshiny)>cmax).and.(cmax>0.)) then
                good(ix,iy) = .false.
                nsat = nsat+1
              endif
              cardd = cards(ix,iy,icshiny)
              card1 = cardd
              if(cardd>0.) then
                card2 = cards(ix,iy,icsky)/cardd
              else
                card2 = 0.
              endif
              cardstd = sqrt(abs(2.*cardd/eperadu)+readvar)
            elseif((modecard==MBLKSKY).or.(modecard==MBSBS) &
                .or.(modecard==MOBJ).or.(modecard==MBLKOBJ) &
                .or.(modecard==MSKY)) then
              if((cards(ix,iy,icblk)>cmax).and.(cmax>0.)) then
                good(ix,iy) = .false.
                nsat = nsat+1
              endif
              if((modecard==MOBJ).or.(modecard==MSKY)) then
                cardd = cards(ix,iy,icsky)
              else
                cardd = cards(ix,iy,icblk)-skyscl*cards(ix,iy,icsky)
              endif
              if(icshiny/=icblk) then
                card1 = cards(ix,iy,icblk)-cards(ix,iy,icshiny)
              else
                card1 = cards(ix,iy,icblk)
              endif
              if(cards(ix,iy,icblk)>0.) then
                card2 = 1.-skyscl*cards(ix,iy,icsky)/cards(ix,iy,icblk)
              else
                card2 = 0.
              endif
              cardstd = sqrt(abs(2.*cards(ix,iy,icblk)/eperadu)+readvar)
            elseif((modecard==MBLKSHINY).and.(icshiny/=icblk)) then
              if((cards(ix,iy,icblk)>cmax).and.(cmax>0.)) then
                good(ix,iy) = .false.
                nsat = nsat+1
              endif
              cardd = cards(ix,iy,icblk)-cards(ix,iy,icshiny)
              card1 = cardd
              if(cards(ix,iy,icblk)>0.) then
                card2 = (cards(ix,iy,icblk)-skyscl*cards(ix,iy,icsky)) &
                  /cards(ix,iy,icblk)
              else
                card2 = 0.
              endif
              cardstd = sqrt(abs(2.*cards(ix,iy,icblk)/eperadu)+readvar)
            else
              if(icshiny/=icblk) then
                print '("*** Unrecognizable cardmode")'
                print '(" Hit RETURN to continue",$)'
                read '(a8)', yn
                write(iupl,'("!!! Unrecognizable cardmode")')
              else
                print '("*** Unusable cardmode without shiny")'
                print '(" Hit RETURN to continue",$)'
                read '(a8)', yn
                write(iupl,'("!!! Unusable cardmode without shiny")')
              endif
              ierr = 9
              if(modeobs==MCAMERA) modecard = modeold
              return
            endif
            if(modecard==MBLK) then
              card1 = cards(ix,iy,icblk)
            elseif(modecard==MSKY) then
              card1 = cards(ix,iy,icsky)
            elseif(modecard==MSHINY) then
              card1 = cards(ix,iy,icshiny)
            endif
            if(card2/=card2) card2 = 0.
            diff(ix,iy) = cardd
            std(ix,iy) = cardstd
            cards(ix,iy,1) = card1
            cards(ix,iy,2) = card2
            cards(ix,iy,3) = cardd
          enddo   
        enddo   
        nxny = nx*ny
        rblsh = (nxny*sumblsh-sumbl*sumsh) &
          /sqrt((nxny*sumbl2-sumbl**2)*(nxny*sumsh2-sumsh**2))
        rsksk = (nxny*sumsksl-sumsk*sumsl) &
          /sqrt((nxny*sumsk2-sumsk**2)*(nxny*sumsl2-sumsl**2))
        if(verbose.or.(rblsh<(0.9)) &
          .or.(rsksk<(0.95)).or.(rsksk>(1.05))) &
          print '(" black-shiny, sky-sky correlation coefficients:", &
          2f8.4)', rblsh,rsksk
        if((rblsh<(0.8)).or.(rsksk<(0.9)).or.(rsksk>(1.1))) &
          print '("*** Likely non-linear or trashed flats")'
        if(nsat>4) then
          print '("*** Warning:",i5," pixels saturated in black"/ &
            " Saturation value =",f10.2)', nsat,cmax
          write(iupl,'("!!!",i5," pixels saturated in black")') nsat
          if(doplot.and.(modecard/=MSHINY)) then
            print '(" Try using cardmode = shiny")'
            if(nsat>32) then
              print '(" Hit RETURN to continue",$)'
              read '(a8)', yn
            endif
          endif
        endif

        if(verbose) then
          print '(" Plotting card diff")'
          if(pause.and.ask) then
            call grey(mx,nx,ny,diff,-1.,-1.,iwin1)
          else
            call grey(mx,nx,ny,diff,0.,0.,iwin1)
          endif
          call waitasec(pause,.true.,diff,mx,nx,ny)
        endif

!	clean and testtort black-shiny (or whatever)
!	testtort finds orders and tests tort with a 2-d fft

        if(modetort>MNONE) then
          nz = 2
          call clean(cards,std,nz,ierr)
!	the array to be testtorted must be in cards(,,1)
          call testtort(cards,ierr)
          if(ierr/=0) then
            print '(" Setting flat = 1")'
            call setone(cards,nc,ierr)
            call setone(flat,1,ierr)
            if(modeobs==MCAMERA) modecard = modeold
            return
          endif
        endif

        if(modecard==MNONE) then
          do iy = 1,ny
            do ix = 1,nx
              do ic = 1,nc
                cards(ix,iy,ic) = 1.
              enddo   
              flat(ix,iy) = 1.
              illum(ix,iy) = 1
            enddo   
          enddo   
          return
        endif

!	now process flat field pair
!	clean diff (normally black-sky)

        nz = 1
        call clean(diff,std,nz,ierr)
        if(verbose) then
          print '(" Plotting cleaned card diff")'
          if(pause.and.ask) then
            call grey(mx,nx,ny,diff,-1.,-1.,iwin1)
          else
            call grey(mx,nx,ny,diff,0.,0.,iwin1)
          endif
          call waitasec(pause,.true.,diff,mx,nx,ny)
          call clean(std,std,nz,ierr)
          print '(" Plotting cleaned card std")'
          if(pause.and.ask) then
            call grey(mx,nx,ny,std,-1.,-1.,iwin1)
          else
            call grey(mx,nx,ny,std,0.,0.,iwin1)
          endif
          call waitasec(pause,.true.,std,mx,nx,ny)
        endif

!	set flat = 0 if diff < 0.25*thrfac*mean to prevent huge flat values
!	note: can't just use illuminated pixels since illum is torted

        sumd = 0.
        sumblk = 0.
        do iy = 1,ny
          do ix = 1,nx
            sumd = sumd+diff(ix,iy)
            sumblk = sumblk+cards(ix,iy,1)
          enddo   
        enddo   
        if(sumd<=0.) then
          print '("*** Mean flat diff <= 0;  Setting flat = 1")'
          print '(" Hit RETURN to continue",$)'
          read '(a8)', yn
          write(iupl,'("!!! Mean flat diff <= 0;  Setting flat = 1")')
        endif
        if(sumd/=sumd) then
          print '("*** Warning: NaN in flat diff")'
          print '(" Hit RETURN to continue",$)'
          read '(a8)', yn
          write(iupl,'("!!! Warning: NaN in flat diff")')
          sumd = 0.
        endif
        if(thrfac<(0.2)) then
          dmin = 0.05*sumd/(nx*ny)
        elseif(thrfac>(1.0)) then
          dmin = 0.25*sumd/(nx*ny)
        else
          dmin = 0.25*thrfac*sumd/(nx*ny)
        endif
        blkmin = 0.1*sumblk/(nx*ny)

!	flat = B_nu(T_card)/diff

        bnut = bnu(waveno0,temp,ierr)
        if((modecard==MOBJ).or.(modecard==MSKY)) bnut = 0.1*bnut
        do iy = 1,ny
          do ix = 1,nx
            if(dmin<=0.) then
              flat(ix,iy) = 1.
            elseif(diff(ix,iy)>dmin) then
              flat(ix,iy) = bnut/diff(ix,iy)
            else
              flat(ix,iy) = 0.
            endif
!	move (blk-sky)/blk to cards(,,1) and blk to cards(,,3) for storespec
            cards(ix,iy,3) = cards(ix,iy,1)
            if(cards(ix,iy,1)>blkmin) then
              cards(ix,iy,1) = cards(ix,iy,2)
            else
              cards(ix,iy,1) = 0.
            endif
          enddo   
        enddo   
        bsmed = quartf(cards(1,1,1),mx*my,3,ierr)
        flatmed = quartf(flat,mx*my,3,ierr)
        print '(" blksky, black med, bnu(T),flat med:",f8.3,2f8.1,f8.3)', &
          bsmed,cmed(icblk),bnut,flatmed

        if(modeinst==MCAMERA) modecard = modeold
        return
      end


      subroutine despike(arr,var,varfac,nz,ierr)
!	look for clouds in frame means,
!	find temporal spikes and replace with average of others,
!	and return variance array after despiking
!	returned variance array is now calculated assuming photon noise
!	unless verbose.and.(stare.or.scan)

        use ius
        use dims
        use modes

        real arr(mx,my,1),var(mx,my)
        real snr(mx,my)
        real sky(mp)
        real nodpa,lores,kmirror,krot
        integer nspike(mp)
        logical good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical scan,shotvar
        character(8) yn
        character(60) comment

        common /nn/ nx,ny,nc,norder,ns,nt
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /iwin/ iwin1,iwin2

!	look at temporal variations in sky to find clouds

        frgain = frtime*abs(gain)
        sixteen = 16./frgain**2
        varmin = 256./(frtime*beamtime*gain**2)
        eperadu = 100.*beamtime
        readvar = (30./eperadu)**2
        scan = (modeobs==MSCAN).or.(modeobs==MMAP)
        shotvar = .not.(verbose.and.(scan.or.(modeobs==MSTARE)))

        if(nz<=2) then
          do iy = 1,ny
            do ix = 1,nx
              var(ix,iy) = readvar
            enddo
          enddo
          return
        endif

  300   sumsky = 0.
        sumsq = 0.
        xnz = 0.
        skymin = 1.0e+31
        do iz = 1,nz
          if((.not.scan).and.(wt(iz)==0.)) cycle
          sum = 0.
          xnxy = 0.
          do iy = 1,ny
            do ix = 1,nx
              if(good(ix,iy)) then
                sum = sum+arr(ix,iy,iz)
                xnxy = xnxy+1.
              endif
            enddo   
          enddo   
          if(xnxy==0.) then
            print '("*** All pixels are bad.  Give up.")'
            print '(" Hit RETURN to continue",$)'
            read '(a8)', yn
            ierr = 9
            return
          endif
          sum = sum/xnxy
          sky(iz) = sum
          sumsky = sumsky+sum
          sumsq = sumsq+sum**2
          if(sum<skymin) skymin = sum
          if(sum/=0.) xnz = xnz+1.
        enddo   
        if(xnz==0.) then
          print '("*** All weights or all pixels = 0 in despike.  Give up.")'
          print '(" Hit RETURN to continue",$)'
          read '(a8)', yn
          ierr = 9
          return
        endif
        avgsky = sumsky/xnz

!	calculate sky noise and check for a trashed frame
        if(spikefac==0.) goto 500
        sumsky = 0.
        sumsq = 0.
        xnz = 0.
        do iz = 1,nz
          if((.not.scan).and.(wt(iz)==0.)) cycle
          skyi = sky(iz)-avgsky
          sumsky = sumsky+skyi
          sumsq = sumsq+skyi**2
          xnz = xnz+1.
        enddo   
        xnz1 = xnz-1.
        xnz = 0.
        sumvar = 0.
        ntrashed = 0
        do iz = 1,nz
          if(scan) then
            skyi = sky(iz)-avgsky
            avg1 = (sumsky-skyi)/xnz1
            var1 = (sumsq-skyi**2)/xnz1-avg1**2
            if((skyi-avg1)**2>(100.*var1)) then
              ntrashed = ntrashed+1
              if(ntrashed==1) then
                print '("*** Fixing trashed frame ",i3)', iz
                if(iz==1) then
                  sky(iz) = sky(iz+1)
                  do iy = 1,ny
                    do ix = 1,nx
                      arr(ix,iy,iz) = arr(ix,iy,iz+1)
                    enddo
                  enddo
                elseif(iz==nz) then
                  sky(iz) = sky(iz-1)
                  do iy = 1,ny
                    do ix = 1,nx
                      arr(ix,iy,iz) = arr(ix,iy,iz-1)
                    enddo
                  enddo
                else
                  sky(iz) = (sky(iz-1)+sky(iz+1))/2.
                  do iy = 1,ny
                    do ix = 1,nx
                      arr(ix,iy,iz) = (arr(ix,iy,iz-1)+arr(ix,iy,iz+1))/2.
                    enddo
                  enddo
                endif
              else
                print '("*** Frame",i3," is trashed.  Skip scan.")', iz
                print '(" Hit RETURN to continue",$)'
                read '(a8)', yn
                ierr = 9
                return
              endif
            elseif((skyi-avg1)**2<(10.*var1)) then
              xnz = xnz+1.
              sumvar = sumvar+var1
            endif
          else
            if(wt(iz)==0.) cycle
            skyi = sky(iz)-avgsky
            avg1 = (sumsky-skyi)/xnz1
            var1 = (sumsq-skyi**2)/xnz1-avg1**2
            if((skyi-avg1)**2>(100.*var1)) then
              print '(" Frame",i3," is trashed.  Skip it.")', iz
              wt(iz) = 0.
              goto 300
            endif
            xnz = xnz+1.
            sumvar = sumvar+var1
          endif
        enddo   

  500   rmssky = sqrt(sumvar/xnz)
        print '(" Sky min, avg, rms: ",3f10.3)', skymin,avgsky,rmssky
        if(dofits) then
          if(skymin>0.) rmssky = rmssky/skymin
          comment = 'sky rms / min'
          call fithreal('SKYNOISE',rmssky,comment,iufith)
!	iuredh may not be open
          write(iuredh,'("skynoise= ",f8.5)') rmssky
        endif
        if(scan.or.(cloud<=0.)) goto 200
        if(cloud<1.) then
          skymin = (1.0+cloud)*skymin
        else
          skymin = cloud
        endif
        ncloud = 0
        do iz = 1,nz
          if((wt(iz)/=0.).and.(sky(iz)>skymin)) then
            wt(iz) = 0.
            ncloud = ncloud+1
          endif
        enddo   
        if(ncloud>0) then
          print '("*** Rejecting ",i2," frames with clouds")', ncloud
          write(iupl,'(" Rejecting ",i2," frames with clouds")') ncloud
          print '(" threshold, cloud =",2es9.2)', skymin,cloud
          print '(8es9.2)', (sky(iz),iz=1,nz)
          xnz = xnz-ncloud
        endif

  200   if(xnz<4.) then
          if(.not.scan) &
            print '(" Can''t despike an array with nz < 4")'
          ierr = -1
          goto 600
        endif

!	look for saturated pixels, calculate variance,
!	then only look for spikes in noisier than average pixels

        amax = satval/frgain
        nsat = 0
        sumavg = 0.
        sumvar = 0.
        do iy = 1,ny
          do ix = 1,nx
            if(.not.good(ix,iy)) then
              var(ix,iy) = 0.
              cycle
            endif
            sum = 0.
            xnz = 0.
            do iz = 1,nz
              if((sky(iz)/=0.).and.(scan.or.(wt(iz)/=0.))) then
                value = arr(ix,iy,iz)*avgsky/sky(iz)
                sum = sum+value
                xnz = xnz+1.
              endif
            enddo   
            avg = sum/xnz
            if((amax>0.).and.(avg>amax)) then
              good(ix,iy) = .false.
              var(ix,iy) = 0.
              nsat = nsat+1
              cycle
            endif
            sumsq = 0.
            sum = 0.
            do iz = 1,nz
              if((sky(iz)/=0.).and.(scan.or.(wt(iz)/=0.))) then
                value = arr(ix,iy,iz)*avgsky/sky(iz)-avg
                sum = sum+value
                sumsq = sumsq+value**2
              endif
            enddo   
!	varfac = 2 if nodding off slit
            vari = varfac*(sumsq/xnz-(sum/xnz)**2)
            if(.not.(vari>=0.)) then
              print '("*** Variance < 0. ",2i4,6es9.2)', &
                ix,iy,sum,sumsq,xnz,vari,avg,avgsky
              print '(" Hit RETURN to continue",$)'
              read '(a8)', yn
            endif
            sumavg = sumavg+avg
            sumvar = sumvar+vari
            var(ix,iy) = vari+sixteen
          enddo   
        enddo
        avgavg = sumavg/xnxy
        avgvar = sumvar/xnxy
        calcvar = readvar+2.*avgavg/eperadu
        print '(" Mean measured, calculated stddev:",2f10.2)', &
          sqrt(avgvar),sqrt(calcvar)
        if(avgvar<varmin) then
          print '(" ** Warning: may be read noise limited.")'
          print '(" RMS std dev, mean / frame =",f6.2,f8.0)', &
            sqrt(avgvar*frtime*beamtime)*abs(gain),avgavg*frgain
          write(iupl, &
            '(" ** Warning: may be read noise limited.")')
          write(iupl,'(" RMS std dev, mean / frame =",f6.2,f8.0)') &
            sqrt(avgvar*frtime*beamtime)*abs(gain),avgavg*frgain
        endif
        avgvar = avgvar+sixteen
        if(nsat>0) print '(i5," saturated pixels found")', nsat

!	a spike is if value-mean' > spikefac*std' (' is w/o value)
        
        nskip = 0
        do iz = 1,nz
          nspike(iz) = 0
          if((.not.scan).and.(wt(iz)==0.)) nskip = nskip+1
        enddo   

        ierr = 0
        if(spikefac<=0.) goto 600
        spsq = spikefac**2
        do iy = 1,ny
          do ix = 1,nx
            if((var(ix,iy)<avgvar).or.(.not.good(ix,iy))) cycle
            sum = 0.
            xnz = 0.
            do iz = 1,nz
              if(scan.or.(wt(iz)/=0.)) then
                value = arr(ix,iy,iz)*avgsky/sky(iz)
                sum = sum+value
                xnz = xnz+1.
              endif
            enddo   
            avg = sum/xnz
            xnz1 = xnz-1.0
            sumv = 0.0
            sumsq = 0.0
            do iz = 1,nz
              if(scan.or.(wt(iz)/=0.)) then
                value = arr(ix,iy,iz)*avgsky/sky(iz)-avg
                sumv = sumv+value
                sumsq = sumsq+value**2
              endif
            enddo   
            avgv = sumv/xnz
            do iz = 1,nz
              if(scan.or.(wt(iz)/=0.)) then
                value = arr(ix,iy,iz)*avgsky/sky(iz)-avg
                avg1 = (sumv-value)/xnz1
                avgsq1 = (sumsq-value**2)/xnz1
!	replace spikes with average of other frames and correct var
                if(((value-avg1)**2)>(spsq*(avgsq1-avg1**2))) then
                  arr(ix,iy,iz) = (avg1+avg)*sky(iz)/avgsky
                  sumv = sumv-value
                  sumsq = sumsq-value**2
                  xnz1 = xnz1-1.
                  var(ix,iy) = varfac*(avgsq1-avg1**2)+sixteen
                  ierr = ierr+1
                  nspike(iz) = nspike(iz)+1
                endif
              endif
            enddo   
          enddo   
        enddo   

        do iz = 1,nz
          if((scan.or.(wt(iz)/=0.)).and.(nspike(iz)>2)) then
            print '(i6," spikes in pair",i3)', nspike(iz),iz
            write(iupl,'(i6," spikes in pair",i3)') nspike(iz),iz
            if(scan.or.(nspike(iz)<=16)) cycle
            nskip = nskip+1
            if((nz-nskip)>=4) then
              print '("*** Skipping pair",i3)', iz
              write(iupl,'(" Skipping pair",i3)') iz
              wt(iz) = 0.
              goto 300
            else
              print '(" *** But skipping it would make nz < 4")'
            endif
          endif
        enddo   

        print '(" Found ",i5," spikes")', ierr

  600   continue
        do iy = 1,ny
        do ix = 1,nx
          snr(ix,iy) = sqrt(var(ix,iy))
        enddo   
        enddo   
        if(verbose) then
          print '(" Plotting measured raw std")'
          call grey(mx,nx,ny,snr,0.,0.,iwin1)
          call waitasec(pause,.true.,snr,mx,nx,ny)
        endif
!	Replace var with calculated shot noise
        do iy = 1,ny
        do ix = 1,nx
          vari = var(ix,iy)
          arri = arr(ix,iy,nz)
          if((vari>0.).and.(arri>0.)) then
            snr(ix,iy) = arri/sqrt(max(readvar,vari))
            if(shotvar) var(ix,iy) = readvar+2.*arri/eperadu
          else
            snr(ix,iy) = 0.
            var(ix,iy) = 0.
          endif
        enddo   
        enddo   
        do iy = 2,ny-1
        do ix = 2,nx-1
           sumsnr = 0.
           do jy = iy-1,iy+1
           do jx = ix-1,ix+1
             sumsnr = sumsnr+snr(jx,jy)
           enddo   
           enddo   
           if(0.*sumsnr==0.) then
             snr(ix,iy) = sumsnr/9.
           else
             snr(ix,iy) = 0.
           endif
        enddo   
        enddo   
        if(verbose) then
          print '(" Plotting SNR image")'
!	  print '(" Plotting measured/calculated std")'
          call grey(mx,nx,ny,snr,0.,0.,iwin1)
          call waitasec(pause,.true.,snr,mx,nx,ny)
        endif

        return

      end


      subroutine debounce(arr,brr,flat,nz,ierr)
!	find the best bounce parameter for each pair to minimize rms difference
!	then debounce by smoothing or shifting arr or brr

        use ius
        use dims
        use modes

        real arr(mx,my,1),brr(mx,my,1)
        real flat(mx,my)
        real abrr(mx,my),diff(mx,my)
        real ai(mp),bi(mp)
        real krot
        logical ok(mx,my),good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical scan

        common /nn/ nx,ny,nc,norder,ns,nt
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder

        if(bounce==0.) return

        print '(" Debouncing data")'
        scan = (modeobs==MSCAN).or.(modeobs==MMAP)
        nzsc = nz-nsky-1

      if(crossdisp) then
!	x-d mode: debounce in x

        nzero = 0
        do iz = 1,nz  !  ends at 100
          ai(iz) = 0.
          bi(iz) = 0.
          if((.not.scan).and.(wt(iz)==0.)) cycle

          do iy = 1,ny
            do ix = 1,nx
              if((ix==1).or.(ix==nx).or.(iy==1).or.(iy==ny))then
                ok(ix,iy) = .false.
              elseif((arr(ix,iy,iz)/=0.).and.(brr(ix,iy,iz)/=0.) &
                .and.good(ix-1,iy).and.good(ix,iy).and.good(ix+1,iy) &
                .and.(flat(ix-1,iy)/=0.).and.(flat(ix,iy)/=0.) &
                .and.(flat(ix+1,iy)/=0.)) then
                ok(ix,iy) = .true.
              else
                ok(ix,iy) = .false.
              endif
            enddo   
          enddo   

!	first fix first derivative bounce by shifting brr
!	calculate variance of arr-brr for three bounce parameters

          ab = abs(bounce)

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+0.25*ab* &
                  (arr(ix+1,iy,iz)-arr(ix-1,iy,iz))
                brri = brr(ix,iy,iz)+0.25*ab* &
                  (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(ix)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          v = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)
                brri = brr(ix,iy,iz)-0.5*ab* &
                  (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(ix)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          u = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)
                brri = brr(ix,iy,iz)+0.5*ab* &
                  (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(ix)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          w = var

!	calculate bounce for minimum variance

          a = u-2.*v+w
          b = u-w
          if(a>(0.0)) then
            ba = 0.5*ab*b/a
          else
            print '("*** Can''t find best 1st derivative bounce", &
              /4x,i8,5es10.3)', iz,u,v,w,ab,0.5*b/a
            if(.not.scan) then
              wt(iz) = 0.
              ai(iz) = -1.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              ai(iz) = 0.
              nzero = nzero+1
            endif
            cycle
          endif
          if(abs(ba)<(ab/100.)) then
            ai(iz) = 0.
            goto 200
          elseif(abs(ba)>(2.*ab)) then
            if(.not.scan) then
              ai(iz) = 0.
              wt(iz) = 0.
              nzero = nzero+1
              goto 100
            elseif(iz<=nzsc) then
              ai(iz) = 0.
              nzero = nzero+1
              goto 100
            else
              ba = 2.*sign(ab,ba)
            endif
          endif
!	skip recalculations for scans
          if(scan) then
            b = 0.
            goto 190
          endif

!	recalculate variance of arr-brr for three bounce parameters

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+0.25*ba* &
                  (arr(ix+1,iy,iz)-arr(ix-1,iy,iz))
                brri = brr(ix,iy,iz)+0.25*ba* &
                  (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(ix)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          u = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)
                brri = brr(ix,iy,iz)+0.5*ba* &
                  (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(ix)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          v = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)
                brri = brr(ix,iy,iz)+ba* &
                  (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(ix)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          w = var

!	calculate bounce for minimum variance and shift brr

          a = u-2.*v+w
          b = u-w
  190     if((a>(0.)).and.(abs(b/a)<(2.))) then
            ba = ba*(1.+0.5*b/a)
          elseif(abs(ba)<(0.1*ab)) then
            continue
          elseif(a>(0.)) then
            if(b>(0.)) then
              ba = 2.*ba
            else
              ba = 0.
            endif
            if(.not.scan) then
              ai(iz) = 0.
              wt(iz) = 0.
              nzero = nzero+1
              goto 100
            elseif(iz<=nzsc) then
              ai(iz) = 0.
              nzero = nzero+1
              goto 100
            endif
          else
            print '("*** Error recalculating 1st derivative bounce", &
              /4x,i8,5es10.3)', iz,u,v,w,ba,0.5*b/a
            if(.not.scan) then
              wt(iz) = 0.
              ai(iz) = -2.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              ai(iz) = 0.
              nzero = nzero+1
            endif
            goto 100
          endif
          ai(iz) = ba
          do iy = 2,ny-1
            do ix = 2,nx-1
              if((brr(ix,iy,iz)/=0.).and.good(ix-1,iy) &
                .and.good(ix,iy).and.good(ix+1,iy)) then
                  diff(ix,iy) = brr(ix,iy,iz)+0.5*ba* &
                    (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
              else
                diff(ix,iy) = brr(ix,iy,iz)
              endif
            enddo   
            do ix = 2,nx-1
              brr(ix,iy,iz) = diff(ix,iy)
            enddo   
          enddo   

!	next fix 2nd derivative bounce by smoothing arr or brr
!	calculate variance of arr-brr for three bounce parameters

  200     if(bounce>0.) goto 100
          ab = abs(bounce)

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+0.5*ab* &
                  (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                brri = brr(ix,iy,iz)+0.5*ab* &
                  (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(ix)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          v = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+ab* &
                  (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                brri = brr(ix,iy,iz)
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(ix)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          u = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)
                brri = brr(ix,iy,iz)+ab* &
                  (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(ix)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          w = var

!	calculate bounce for minimum variance

          a = u-2.*v+w
          b = u-w
          if(a>(0.0)) then
            ba = 0.5*ab*b/a
          else
            print '("*** Can''t find best 2nd derivative bounce", &
              /4x,i8,5es10.3)', iz,u,v,w,ab,0.5*b/a
            if(.not.scan) then
              wt(iz) = 0.
              bi(iz) = -1.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              bi(iz) = 0.
              nzero = nzero+1
            endif
            goto 100
          endif
          if(abs(ba)<(ab/100.)) then
            bi(iz) = 0.
            goto 100
          elseif(abs(ba)>(4.*ab)) then
            if(.not.scan) then
              bi(iz) = 0.
              wt(iz) = 0.
              nzero = nzero+1
              goto 100
            elseif(iz<=nzsc) then
              bi(iz) = 0.
              nzero = nzero+1
              goto 100
            else
              ba = 4.*sign(ab,ba)
            endif
          endif
!	skip recalculations for scans
          if(scan) then
            b = 0.
            goto 290
          endif

!	recalculate variance of arr-brr for three bounce parameters

          ab = abs(ba)

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+ab* &
                  (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                brri = brr(ix,iy,iz)+ab* &
                  (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(ix)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          u = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                if(ba<0.) then
                  arri = arr(ix,iy,iz)+1.5*ab* &
                   (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                  brri = brr(ix,iy,iz)+0.5*ab* &
                   (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                elseif(ba>0.) then
                  arri = arr(ix,iy,iz)+0.5*ab* &
                   (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                  brri = brr(ix,iy,iz)+1.5*ab* &
                   (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                endif
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(ix)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          v = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                if(ba<0.) then
                  arri = arr(ix,iy,iz)+2.*ab* &
                   (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                  brri = brr(ix,iy,iz)
                elseif(ba>0.) then
                  arri = arr(ix,iy,iz)
                  brri = brr(ix,iy,iz)+2.*ab* &
                   (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                endif
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(ix)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          w = var

!	calculate bounce for minimum variance and debounce arr or brr

          a = u-2.*v+w
          b = u-w
  290     if((a>0.).and.(abs(b/a)<2.)) then
            ba = ba*(1.+0.5*b/a)
            bi(iz) = ba
          elseif(abs(ba)<(0.1*abs(bounce))) then
            bi(iz) = ba
          else
            print '("*** Error recalculating 2nd derivative bounce", &
              /4x,i8,5es10.3)', iz,u,v,w,ba,0.5*b/a
            if(.not.scan) then
              wt(iz) = 0.
              bi(iz) = -2.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              bi(iz) = 0.
              nzero = nzero+1
            endif
            goto 100
          endif
          ab = abs(ba)
          do iy = 2,ny-1
            do ix = 2,nx-1
              if(ba<0.) then
                abrr(ix,iy) = arr(ix,iy,iz)+ab* &
                 (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
              elseif(ba>0.) then
                abrr(ix,iy) = brr(ix,iy,iz)+ab* &
                 (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
              endif
            enddo   
          enddo   
          do iy = 2,ny-1
            do ix = 2,nx-1
              if(ok(ix,iy)) then
                if(ba<0.) then
                  arr(ix,iy,iz) = abrr(ix,iy)
                elseif(ba>0.) then
                  brr(ix,iy,iz) = abrr(ix,iy)
                endif
              endif
            enddo   
          enddo   

  100   enddo   

      else
!	long-slit mode: debounce in y (deleted)
	goto 301

        nzero = 0
        do iz = 1,nz  !  ends at 300
          ai(iz) = 0.
          bi(iz) = 0.
          if((.not.scan).and.(wt(iz)==0.)) cycle

          do iy = 1,ny
            do ix = 1,nx
              if((ix==1).or.(ix==nx).or.(iy==1).or.(iy==ny))then
                ok(ix,iy) = .false.
              elseif((arr(ix,iy,iz)/=0.).and.(brr(ix,iy,iz)/=0.) &
                .and.good(ix,iy-1).and.good(ix,iy).and.good(ix,iy+1) &
                .and.(flat(ix,iy-1)/=0.).and.(flat(ix,iy)/=0.) &
                .and.(flat(ix,iy+1)/=0.)) then
                ok(ix,iy) = .true.
              else
                ok(ix,iy) = .false.
              endif
            enddo   
          enddo   

!	first fix first derivative bounce by shifting brr
!	calculate variance of arr-brr for three bounce parameters

          ab = abs(bounce)

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+0.25*ab* &
                  (arr(ix,iy+1,iz)-arr(ix,iy-1,iz))
                brri = brr(ix,iy,iz)+0.25*ab* &
                  (brr(ix,iy+1,iz)-brr(ix,iy-1,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do iy = 1,ny
            if(illx(iy)/=1) cycle
            sum = 0.
            do ix = 1,nx
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          v = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)
                brri = brr(ix,iy,iz)-0.5*ab* &
                  (brr(ix,iy+1,iz)-brr(ix,iy-1,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do iy = 1,ny
            if(illx(iy)/=1) cycle
            sum = 0.
            do ix = 1,nx
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          u = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)
                brri = brr(ix,iy,iz)+0.5*ab* &
                  (brr(ix,iy+1,iz)-brr(ix,iy-1,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do iy = 1,ny
            if(illx(iy)/=1) cycle
            sum = 0.
            do ix = 1,nx
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          w = var

!	calculate bounce for minimum variance

          a = u-2.*v+w
          b = u-w
          if(a>(0.0)) then
            ba = 0.5*ab*b/a
          else
            print '("*** Can''t find best 1st derivative bounce", &
              /4x,i8,5es10.3)', iz,u,v,w,ab,0.5*b/a
            if(.not.scan) then
              wt(iz) = 0.
              ai(iz) = -1.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              ai(iz) = 0.
              nzero = nzero+1
            endif
            goto 300
          endif
          if(abs(ba)<(ab/100.)) then
            ai(iz) = 0.
            goto 400
          elseif(abs(ba)>(2.*ab)) then
            if(.not.scan) then
              wt(iz) = 0.
              ai(iz) = -1.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              ai(iz) = 0.
              nzero = nzero+1
            endif
            goto 300
!	skip recalculations for scans
          elseif(scan) then
            b = 0.
            goto 390
          endif

!	recalculate variance of arr-brr for three bounce parameters

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+0.25*ba* &
                  (arr(ix,iy+1,iz)-arr(ix,iy-1,iz))
                brri = brr(ix,iy,iz)+0.25*ba* &
                  (brr(ix,iy+1,iz)-brr(ix,iy-1,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do iy = 1,ny
            if(illx(iy)/=1) cycle
            sum = 0.
            do ix = 1,nx
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          u = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)
                brri = brr(ix,iy,iz)+0.5*ba* &
                  (brr(ix,iy+1,iz)-brr(ix,iy-1,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do iy = 1,ny
            if(illx(iy)/=1) cycle
            sum = 0.
            do ix = 1,nx
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          v = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)
                brri = brr(ix,iy,iz)+ba* &
                  (brr(ix,iy+1,iz)-brr(ix,iy-1,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do iy = 1,ny
            if(illx(iy)/=1) cycle
            sum = 0.
            do ix = 1,nx
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          w = var

!	calculate bounce for minimum variance and shift brr

          a = u-2.*v+w
          b = u-w
  390     if((a>(0.)).and.(abs(b/a)<(2.))) then
            ba = ba*(1.+0.5*b/a)
          elseif(abs(ba)>(0.1*ab)) then
            print '("*** Error recalculating 1st derivative bounce", &
              /4x,i8,5es10.3)', iz,u,v,w,ba,0.5*b/a
            if(.not.scan) then
              wt(iz) = 0.
              ai(iz) = -2.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              ai(iz) = 0.
              nzero = nzero+1
            endif
            goto 300
          endif
          ai(iz) = ba
          do iy = 2,ny-1
            do ix = 2,nx-1
              if((brr(ix,iy,iz)/=0.).and.good(ix,iy-1) &
                .and.good(ix,iy).and.good(ix,iy+1)) then
                  diff(ix,iy) = brr(ix,iy,iz)+0.5*ba* &
                    (brr(ix,iy+1,iz)-brr(ix,iy-1,iz))
              else
                diff(ix,iy) = brr(ix,iy,iz)
              endif
            enddo   
          enddo   
          do iy = 2,ny-1
            do ix = 2,nx-1
              if((brr(ix,iy,iz)/=0.).and.good(ix,iy-1) &
                .and.good(ix,iy).and.good(ix,iy+1)) &
                  brr(ix,iy,iz) = diff(ix,iy)
            enddo   
          enddo   

!	next fix 2nd derivative bounce by smoothing arr or brr
!	calculate variance of arr-brr for three bounce parameters

  400     if(bounce>0.) goto 300
          ab = abs(bounce)

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+0.5*ab* &
                  (arr(ix,iy-1,iz)-2.*arr(ix,iy,iz)+arr(ix,iy+1,iz))
                brri = brr(ix,iy,iz)+0.5*ab* &
                  (brr(ix,iy-1,iz)-2.*brr(ix,iy,iz)+brr(ix,iy+1,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do iy = 1,ny
            if(illx(iy)/=1) cycle
            sum = 0.
            do ix = 1,nx
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          v = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+ab* &
                  (arr(ix,iy-1,iz)-2.*arr(ix,iy,iz)+arr(ix,iy+1,iz))
                brri = brr(ix,iy,iz)
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do iy = 1,ny
            if(illx(iy)/=1) cycle
            sum = 0.
            do ix = 1,nx
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          u = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)
                brri = brr(ix,iy,iz)+ab* &
                  (brr(ix,iy-1,iz)-2.*brr(ix,iy,iz)+brr(ix,iy+1,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do iy = 1,ny
            if(illx(iy)/=1) cycle
            sum = 0.
            do ix = 1,nx
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          w = var

!	calculate bounce for minimum variance

          a = u-2.*v+w
          b = u-w
          if(a>(0.0)) then
            ba = 0.5*ab*b/a
          else
            print '("*** Can''t find best 2nd derivative bounce", &
              /4x,i8,5es10.3)', iz,u,v,w,ab,0.5*b/a
            if(.not.scan) then
              wt(iz) = 0.
              bi(iz) = -1.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              bi(iz) = 0.
              nzero = nzero+1
            endif
            goto 300
          endif
          if(abs(ba)<(ab/100.)) then
            bi(iz) = 0.
            goto 300
          elseif(abs(ba)>(4.*ab)) then
            if(.not.scan) then
              wt(iz) = 0.
              bi(iz) = -1.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              bi(iz) = 0.
              nzero = nzero+1
            endif
            goto 300
          elseif(scan) then
!	skip recalculations for scans
            b = 0.
            goto 490
          endif

!	recalculate variance of arr-brr for three bounce parameters

          ab = abs(ba)

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+ab* &
                  (arr(ix,iy-1,iz)-2.*arr(ix,iy,iz)+arr(ix,iy+1,iz))
                brri = brr(ix,iy,iz)+ab* &
                  (brr(ix,iy-1,iz)-2.*brr(ix,iy,iz)+brr(ix,iy+1,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do iy = 1,ny
            if(illx(iy)/=1) cycle
            sum = 0.
            do ix = 1,nx
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          u = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                if(ba<0.) then
                  arri = arr(ix,iy,iz)+1.5*ab* &
                   (arr(ix,iy-1,iz)-2.*arr(ix,iy,iz)+arr(ix,iy+1,iz))
                  brri = brr(ix,iy,iz)+0.5*ab* &
                   (brr(ix,iy-1,iz)-2.*brr(ix,iy,iz)+brr(ix,iy+1,iz))
                elseif(ba>0.) then
                  arri = arr(ix,iy,iz)+0.5*ab* &
                   (arr(ix,iy-1,iz)-2.*arr(ix,iy,iz)+arr(ix,iy+1,iz))
                  brri = brr(ix,iy,iz)+1.5*ab* &
                   (brr(ix,iy-1,iz)-2.*brr(ix,iy,iz)+brr(ix,iy+1,iz))
                endif
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do iy = 1,ny
            if(illx(iy)/=1) cycle
            sum = 0.
            do ix = 1,nx
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          v = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                if(ba<0.) then
                  arri = arr(ix,iy,iz)+2.*ab* &
                   (arr(ix,iy-1,iz)-2.*arr(ix,iy,iz)+arr(ix,iy+1,iz))
                  brri = brr(ix,iy,iz)
                elseif(ba>0.) then
                  arri = arr(ix,iy,iz)
                  brri = brr(ix,iy,iz)+2.*ab* &
                   (brr(ix,iy-1,iz)-2.*brr(ix,iy,iz)+brr(ix,iy+1,iz))
                endif
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do iy = 1,ny
            if(illx(iy)/=1) cycle
            sum = 0.
            do ix = 1,nx
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          w = var

!	calculate bounce for minimum variance and debounce arr or brr

          a = u-2.*v+w
          b = u-w
  490     if((a>0.).and.(abs(b/a)<2.)) then
            ba = ba*(1.+0.5*b/a)
            bi(iz) = ba
            ab = abs(ba)
            do iy = 2,ny-1
              do ix = 2,nx-1
                if(ba<0.) then
                  abrr(ix,iy) = arr(ix,iy,iz)+ab* &
                   (arr(ix,iy-1,iz)-2.*arr(ix,iy,iz)+arr(ix,iy+1,iz))
                elseif(ba>0.) then
                  abrr(ix,iy) = brr(ix,iy,iz)+ab* &
                   (brr(ix,iy-1,iz)-2.*brr(ix,iy,iz)+brr(ix,iy+1,iz))
                endif
              enddo   
            enddo   
            do iy = 1,ny
              do ix = 1,nx
                if(ok(ix,iy)) then
                  if(ba<0.) then
                    arr(ix,iy,iz) = abrr(ix,iy)
                  elseif(ba>0.) then
                    brr(ix,iy,iz) = abrr(ix,iy)
                  endif
                endif
              enddo   
            enddo   
          elseif(abs(ba)<(0.1*abs(bounce))) then
            bi(iz) = ba
            ab = abs(ba)
            do iy = 2,ny-1
              do ix = 2,nx-1
                if(ba<0.) then
                  abrr(ix,iy) = arr(ix,iy,iz)+ab* &
                   (arr(ix,iy-1,iz)-2.*arr(ix,iy,iz)+arr(ix,iy+1,iz))
                elseif(ba>0.) then
                  abrr(ix,iy) = brr(ix,iy,iz)+ab* &
                   (brr(ix,iy-1,iz)-2.*brr(ix,iy,iz)+brr(ix,iy+1,iz))
                endif
              enddo   
            enddo   
            do iy = 1,ny
              do ix = 1,nx
                if(ok(ix,iy)) then
                  if(ba<0.) then
                    arr(ix,iy,iz) = abrr(ix,iy)
                  elseif(ba>0.) then
                    brr(ix,iy,iz) = abrr(ix,iy)
                  endif
                endif
              enddo   
            enddo   
          else
            print '("*** Error recalculating 2nd derivative bounce", &
              /4x,i8,5es10.3)', iz,u,v,w,ba,0.5*b/a
            if(.not.scan) then
              wt(iz) = 0.
              bi(iz) = -2.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              bi(iz) = 0.
              nzero = nzero+1
            endif
          endif

  300   enddo   
  301   continue

!	long-slit mode: debounce in x

        nzero = 0
        do iz = 1,nz  !  ends at 500
          ai(iz) = 0.
          bi(iz) = 0.
          if((.not.scan).and.(wt(iz)==0.)) cycle

          do iy = 1,ny
            do ix = 1,nx
              ok(ix,iy) = (ix>4).and.(ix<(nx-3)) &
                .and.(iy>1).and.(iy<ny) &
                .and.good(ix-1,iy).and.good(ix,iy).and.good(ix+1,iy) &
                .and.(flat(ix-1,iy)/=0.).and.(flat(ix,iy)/=0.) &
                .and.(flat(ix+1,iy)/=0.)
            enddo   
          enddo   

!	first fix first derivative bounce by shifting brr
!	calculate variance of arr-brr for three bounce parameters

          ab = abs(bounce)

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+0.25*ab* &
                  (arr(ix+1,iy,iz)-arr(ix-1,iy,iz))
                brri = brr(ix,iy,iz)+0.25*ab* &
                  (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          v = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+0.25*ab* &
                  (arr(ix+1,iy,iz)-arr(ix-1,iy,iz))
                brri = brr(ix,iy,iz)-0.25*ab* &
                  (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          u = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)-0.25*ab* &
                  (arr(ix+1,iy,iz)-arr(ix-1,iy,iz))
                brri = brr(ix,iy,iz)+0.25*ab* &
                  (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          w = var

!	calculate bounce for minimum variance

          a = u-2.*v+w
          b = u-w
          if(a>(0.0)) then
            ba = 0.5*ab*b/a
          else
            print '("*** Can''t find best 1st derivative bounce", &
              /4x,i8,5es10.3)', iz,u,v,w,ab,0.5*b/a
            if(.not.scan) then
              wt(iz) = 0.
              ai(iz) = -1.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              ai(iz) = 0.
              nzero = nzero+1
            endif
            goto 500
          endif
          if(abs(ba)<(ab/100.)) then
            ai(iz) = 0.
            goto 600
          elseif(abs(ba)>(2.*ab)) then
            if(.not.scan) then
              wt(iz) = 0.
              ai(iz) = -1.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              ai(iz) = 0.
              nzero = nzero+1
            endif
            goto 500
!	skip recalculations for scans
          elseif(scan) then
            b = 0.
            goto 590
          endif

!	recalculate variance of arr-brr for three bounce parameters

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+0.5*ba* &
                  (arr(ix+1,iy,iz)-arr(ix-1,iy,iz))
                brri = brr(ix,iy,iz)+0.5*ba* &
                  (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          u = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+0.25*ba* &
                  (arr(ix+1,iy,iz)-arr(ix-1,iy,iz))
                brri = brr(ix,iy,iz)+0.75*ba* &
                  (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          v = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)-0.5*ba* &
                  (arr(ix+1,iy,iz)-arr(ix-1,iy,iz))
                brri = brr(ix,iy,iz)+0.5*ba* &
                  (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          w = var

!	calculate bounce for minimum variance and shift brr

          a = u-2.*v+w
          b = u-w
  590     if((a>(0.)).and.(abs(b/a)<(2.))) then
            ba = ba*(1.+0.5*b/a)
          elseif(abs(ba)>(0.1*ab)) then
            print '("*** Error recalculating 1st derivative bounce", &
              /4x,i8,5es10.3)', iz,u,v,w,ba,0.5*b/a
            if(.not.scan) then
              wt(iz) = 0.
              ai(iz) = -2.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              ai(iz) = 0.
              nzero = nzero+1
            endif
            goto 500
          endif
          ai(iz) = ba
          var = 0.
          do ix = 2,nx-1
            do iy = 2,ny-1
              if((brr(ix,iy,iz)/=0.).and.good(ix-1,iy) &
                .and.good(ix,iy).and.good(ix+1,iy)) then
                  diff(ix,iy) = brr(ix,iy,iz)+0.5*ba* &
                    (brr(ix+1,iy,iz)-brr(ix-1,iy,iz))
              else
                diff(ix,iy) = brr(ix,iy,iz)
              endif
            enddo   
            sum = 0.
            do iy = 2,ny-1
              if((brr(ix,iy,iz)/=0.).and.good(ix-1,iy) &
                .and.good(ix,iy).and.good(ix+1,iy)) &
                  brr(ix,iy,iz) = diff(ix,iy)
              if(ok(ix,iy)) &
                sum = sum+arr(ix,iy,iz)-brr(ix,iy,iz)
            enddo   
            var = var+sum**2
          enddo   

!	next fix 2nd derivative bounce by smoothing arr or brr
!	calculate variance of arr-brr for three bounce parameters

  600     if(bounce>0.) goto 500
          ab = abs(bounce)

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+0.5*ab* &
                  (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                brri = brr(ix,iy,iz)+0.5*ab* &
                  (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          v = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+ab* &
                  (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                brri = brr(ix,iy,iz)
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          u = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)
                brri = brr(ix,iy,iz)+ab* &
                  (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            if(illx(iy)/=1) cycle
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          w = var

!	calculate bounce for minimum variance

          a = u-2.*v+w
          b = u-w
          if(a>(0.0)) then
            ba = 0.5*ab*b/a
          else
            print '("*** Can''t find best 2nd derivative bounce", &
              /4x,i8,5es10.3)', iz,u,v,w,ab,0.5*b/a
            if(.not.scan) then
              wt(iz) = 0.
              bi(iz) = -1.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              bi(iz) = 0.
              nzero = nzero+1
            endif
            goto 500
          endif
          if(abs(ba)<(ab/100.)) then
            bi(iz) = 0.
            goto 500
          elseif(abs(ba)>(4.*ab)) then
            if(.not.scan) then
              wt(iz) = 0.
              bi(iz) = -1.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              bi(iz) = 0.
              nzero = nzero+1
            endif
            goto 500
          elseif(scan) then
!	skip recalculations for scans
            b = 0.
            goto 690
          endif

!	recalculate variance of arr-brr for three bounce parameters

          ab = abs(ba)

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                arri = arr(ix,iy,iz)+ab* &
                  (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                brri = brr(ix,iy,iz)+ab* &
                  (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          u = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                if(ba<0.) then
                  arri = arr(ix,iy,iz)+1.5*ab* &
                   (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                  brri = brr(ix,iy,iz)+0.5*ab* &
                   (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                elseif(ba>0.) then
                  arri = arr(ix,iy,iz)+0.5*ab* &
                   (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                  brri = brr(ix,iy,iz)+1.5*ab* &
                   (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                endif
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          v = var

          do iy = 1,ny
            do ix = 1,nx
              if(ok(ix,iy)) then
                if(ba<0.) then
                  arri = arr(ix,iy,iz)+2.*ab* &
                   (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                  brri = brr(ix,iy,iz)
                elseif(ba>0.) then
                  arri = arr(ix,iy,iz)
                  brri = brr(ix,iy,iz)+2.*ab* &
                   (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                endif
                diff(ix,iy) = arri-brri
              else
                diff(ix,iy) = 0.
              endif
            enddo   
          enddo   
          call qtort(diff,1,ierr)
          var = 0.
          do ix = 1,nx
            sum = 0.
            do iy = 1,ny
              sum = sum+diff(ix,iy)
            enddo   
            var = var+sum**2
          enddo   
          w = var

!	calculate bounce for minimum variance and debounce arr or brr

          a = u-2.*v+w
          b = u-w
  690     if((a>0.).and.(abs(b/a)<2.)) then
            ba = ba*(1.+0.5*b/a)
            bi(iz) = ba
            ab = abs(ba)
            do iy = 2,ny-1
              do ix = 2,nx-1
                if(ba<0.) then
                  abrr(ix,iy) = arr(ix,iy,iz)+ab* &
                   (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                elseif(ba>0.) then
                  abrr(ix,iy) = brr(ix,iy,iz)+ab* &
                   (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                endif
              enddo   
            enddo   
            do iy = 2,ny-1
              do ix = 2,nx-1
                if(ok(ix,iy)) then
                  if(ba<0.) then
                    arr(ix,iy,iz) = abrr(ix,iy)
                  elseif(ba>0.) then
                    brr(ix,iy,iz) = abrr(ix,iy)
                  endif
                endif
              enddo   
            enddo   
          elseif(abs(ba)<(0.1*abs(bounce))) then
            bi(iz) = ba
            ab = abs(ba)
            do iy = 2,ny-1
              do ix = 2,nx-1
                if(ba<0.) then
                  abrr(ix,iy) = arr(ix,iy,iz)+ab* &
                   (arr(ix-1,iy,iz)-2.*arr(ix,iy,iz)+arr(ix+1,iy,iz))
                elseif(ba>0.) then
                  abrr(ix,iy) = brr(ix,iy,iz)+ab* &
                   (brr(ix-1,iy,iz)-2.*brr(ix,iy,iz)+brr(ix+1,iy,iz))
                endif
              enddo   
            enddo   
            do iy = 2,ny-1
              do ix = 2,nx-1
                if(ok(ix,iy)) then
                  if(ba<0.) then
                    arr(ix,iy,iz) = abrr(ix,iy)
                  elseif(ba>0.) then
                    brr(ix,iy,iz) = abrr(ix,iy)
                  endif
                endif
              enddo   
            enddo   
          else
            print '("*** Error recalculating 2nd derivative bounce", &
              /4x,i8,5es10.3)', iz,u,v,w,ba,0.5*b/a
            if(.not.scan) then
              wt(iz) = 0.
              bi(iz) = -2.
              nzero = nzero+1
            elseif(iz<=nzsc) then
              bi(iz) = 0.
              nzero = nzero+1
            endif
          endif

  500   enddo

      endif

        print '(" First derivative bounce parameters:")'
        print '(8f8.3)', (ai(iz),iz=1,nz)
        write(iupl,'(" First derivative bounce parameters:")')
        write(iupl,'(8f8.3)') (ai(iz),iz=1,nz)
        if(bounce<0.) then
          print '(" Second derivative bounce parameters:")'
          print '(8f8.3)', (bi(iz),iz=1,nz)
          write(iupl,'(" Second derivative bounce parameters:")')
          write(iupl,'(8f8.3)') (bi(iz),iz=1,nz)
        endif
        ierr = 0
        if(nzero>0) then
          if(scan) then
            print '(i3," pairs exceeded bounce limit")', nzero
            ierr = nzero
          else
            print '(" Setting weight = 0 for",i3," pairs")', nzero
          endif
        endif

        return
      end


      subroutine diffarr(a,b,d,avar,bvar,std,nz,ierr)
!	subtract beams and calculate noise

        use ius
        use dims

        real a(mx,my,1),b(mx,my,1),d(mx,my,1)
        real avar(mx,my),bvar(mx,my),std(mx,my)
        logical stare,scan
        real nodpa,lores,kmirror,krot

        common /nn/ nx,ny,nc,norder,ns,nt
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder

        stare = (ierr==1)
        scan = (ierr==2)
        ierr = 0

        dmin = 1./(beamtime*abs(gain))
        nbad = 0
        do iz = 1,nz
          do iy = 1,ny
            do ix = 1,nx
              diff = a(ix,iy,iz)-b(ix,iy,iz)
              if(0.*diff/=0.) then
                if(nbad<4) print '(" In diffarr:",3i8,2es10.2)', &
                  ix,iy,iz,a(ix,iy,iz),b(ix,iy,iz)
                nbad = nbad+1
                d(ix,iy,iz) = 0.
              elseif(diff==0.) then
                d(ix,iy,iz) = dmin
              else
                d(ix,iy,iz) = diff
              endif
            enddo   
          enddo   
        enddo   
        if(nbad>0) print '(" nbad =",i10)', nbad

        dmin = 16./(beamtime*gain)**2
        eperadu = 100.*beamtime
        readvar = (30./eperadu)**2
        izsky = min(3,nz)
        do iy = 1,ny
          do ix = 1,nx
            if(stare) then
              dvar = avar(ix,iy)
            elseif(scan) then
              if(nz<mp) d(ix,iy,nz+1) = b(ix,iy,izsky)
              dvar = 2.*a(ix,iy,izsky)/eperadu
            else
!	      dvar = avar(ix,iy)+bvar(ix,iy)
              dvar = 2.*min(avar(ix,iy),bvar(ix,iy))
            endif
            dvar = max(dvar,dmin)
            std(ix,iy) = sqrt(dvar)
          enddo   
        enddo   

        return
      end


      subroutine calibrate(arr,flat,nz,ierr)
!	multiply array by flat to calibrate

        use ius
        use dims
        use modes

        real arr(mx,my,1),flat(mx,my)
        real nodpa,lores,kmirror,krot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        character(60) comment

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf

        if(modecard==MNONE) then
          if(ierr==0) then
            write(iuredh,'("calibrat= F")')
            comment = 'data were not intensity calibrated'
            if(dofits) call fithlog('CALIBRAT',.false.,comment,iufith)
          endif
          return
        elseif(nz>1) then
          if(ierr==0) then
            write(iuredh,'("calibrat= T")')
            comment = 'data were intensity calibrated with blackbody'
            if(dofits) call fithlog('CALIBRAT',.true.,comment,iufith)
          endif
        endif

        do iz = 1,nz
          do iy = 1,ny
            do ix = 1,nx
              arr(ix,iy,iz) = arr(ix,iy,iz)*flat(ix,iy)
            enddo   
          enddo   
        enddo

        ierr = 0
        return

      end


      subroutine clean(arr,std,nz,ierr)
!	interpolate over bad and noisy pixels

        use ius
        use dims
        use modes

        real arr(mx,my,1),std(mx,my)
        real nodpa,lores,kmirror,krot
        logical good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical hdopen
        character(60) object,feature,obsmode,instmode
        character(60) objtype,piname,pid,note,weather,warning,version
        character(60) comment

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
                piname,pid,note,weather,warning,version
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

!	look for noisy and too quiet pixels

        if((modeobs==MSCAN).or.(stdfac<=0.)) goto 100
        sumstd = 0.
        do iy = 1,ny
          do ix = 1,nx
            if(good(ix,iy)) then
              sumstd = sumstd+std(ix,iy)
            endif
          enddo   
        enddo   
        if(sumstd/=sumstd) print '("*** Warning: NaN in std")'
        if(sumstd<=0.) goto 100
        stdmax = stdfac*sumstd/(nx*ny)
        stdmin = sumstd/(16*nx*ny)
        ibad = 0
        do iy = 1,ny
          do ix = 1,nx
            if(good(ix,iy).and.(std(ix,iy)>stdmax)) then
              good(ix,iy) = .false.
              ibad = ibad+1
            elseif(good(ix,iy).and.(std(ix,iy)<stdmin)) then
              std(ix,iy) = stdmin
            endif
          enddo   
        enddo   

        print '(" Found ",i5," noisy pixels")', ibad
        write(iupl,'(" Found ",i5," noisy pixels")') ibad

  100   continue

!	fix bad pixels by interpolation between neighbors

        ierr = 0
        ibad = 0
        do iy = 1,ny
          do ix = 1,nx
            if(.not.good(ix,iy)) then
              ibad = ibad+1
!	first try y neighbors
              if((iy>1).and.(iy<ny).and. &
                good(ix,iy-1).and.good(ix,iy+1)) then
                do iz = 1,nz
                  arr(ix,iy,iz) = &
                    (arr(ix,iy-1,iz)+arr(ix,iy+1,iz))/2.
                enddo
                std(ix,iy) = &
                  sqrt(std(ix,iy-1)**2+std(ix,iy+1)**2)
!	next try x neighbors
              elseif((ix>1).and.(ix<nx).and. &
                good(ix-1,iy).and.good(ix+1,iy)) then
                do iz = 1,nz
                  arr(ix,iy,iz) = &
                    (arr(ix-1,iy,iz)+arr(ix+1,iy,iz))/2.
                enddo
                std(ix,iy) = &
                  sqrt(std(ix-1,iy)**2+std(ix+1,iy)**2)
              else
!	next look farther in y
                iy1 = iy
  210           iy1 = iy1-1
                if(iy1<1) goto 220
                if(.not.good(ix,iy1)) goto 210
                iy2 = iy
  212           iy2 = iy2+1
                if(iy2>ny) goto 220
                if(.not.good(ix,iy2)) goto 212
!	interpolate in y
                do iz = 1,nz
                  arr(ix,iy,iz) = ((iy2-iy)*arr(ix,iy1,iz) &
                    +(iy-iy1)*arr(ix,iy2,iz))/(iy2-iy1)
                enddo
                std(ix,iy) = &
                  sqrt(std(ix,iy1)**2+std(ix,iy2)**2)
                cycle
!	finally look farther in x
  220           ix1 = ix
  222           ix1 = ix1-1
                if((ix1>0).and.(.not.good(ix1,iy))) goto 222
                ix2 = ix
  224           ix2 = ix2+1
                if((ix2<=nx).and.(.not.good(ix2,iy))) goto 224
!	if on an edge
                if((ix1<1).and.(ix2>nx)) then
                  if(ierr<8) print '(" Can''t clean ",2i4)', ix,iy
                  do iz = 1,nz
                    arr(ix,iy,iz) = 0.
                  enddo
                  std(ix,iy) = 0.
                  ierr = ierr+1
                  if(ierr>nx) goto 200
                elseif(ix1<1) then
                  do iz = 1,nz
                    arr(ix,iy,iz) = arr(ix2,iy,iz)
                  enddo
                  std(ix,iy) = 2.*std(ix2,iy)
                elseif(ix2>nx) then
                  do iz = 1,nz
                    arr(ix,iy,iz) = arr(ix1,iy,iz)
                  enddo
                  std(ix,iy) = 2.*std(ix1,iy)
                else
!	interpolate in x
                  do iz = 1,nz
                    arr(ix,iy,iz) = ((ix-ix2)*arr(ix1,iy,iz) &
                      +(ix-ix1)*arr(ix2,iy,iz))/(ix2-ix1)
                  enddo
                  std(ix,iy) = &
                    sqrt(std(ix1,iy)**2+std(ix2,iy)**2)
                endif
              endif
            endif
          enddo   
        enddo   

  200   inquire(unit=iuredh,opened=hdopen)
        if(ierr>0) then
          print '(i4," could not be cleaned")', ierr
          write(iupl,'(i4," could not be cleaned")') ierr
          if(hdopen.and.(nz>1)) then
            write(iuredh,'("clean   = F")')
            if(dofits) then
              comment = 'could not clean bad pixels'
              call fithlog('CLEAN   ',.false.,comment,iufith)
              if(dosum) call fithlog('CLEAN   ',.false.,comment,iufish)
            endif
          endif
          if(ierr<=nx) ierr = 0
        else
          if(hdopen.and.(nz>1) &
            .and.(modeobs>MFLAT).and.(modeobs<MSCAN)) then
            write(iuredh,'("clean   = T")')
            comment = 'bad pixels cleaned'
            if(dofits) call fithlog('CLEAN   ',.true.,comment,iufith)
          endif
        endif

        return
      end


      subroutine tort(a,nz,ierr)
!	remove optical distortions
!	jun03 version with integral order spacing
!	may05 add sinc weighting option

        use ius
        use dims
        use modes
        use consts

        real a(mx,my,1),b(mx,my)
        real nodpa,lores,kmirror,krot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical good
        logical hdopen
        character(60) comment

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /iwin/ iwin1,iwin2

        if(modetort==MNONE) return

!	calculate mapping of undistorted (x,y) onto arr (u,v)
!	then interpolate from u,v onto x,y

        if(crossdisp) then
          hrslrot = -slitrot+krot
!	include echelon smile and slit rotation here
          hrskew = 2.*hrg*abs(hrr)+tan(hrslrot)
          hrsmile = -abs(hrr)*pixelwd/hrfl
!	should include the variation of hrr along an order.  some day.
!	include xd smile and spectrum rotation by k mirror here
          xdskew = 2.*xdg*xdr+tan(krot)
          xdsmile = -xdr*pixelwd/xdfl
!	xd smile seems to be overcorrected.  put in fudge factor.
          xdsmile = 0.7*xdsmile
          xddisp = (xdr*xdfl)/(abs(hrr)*hrfl)
          xdnlin = -(xdr+1./(2.*xdr))*(pixelwd/xdfl)/2.
!	xdskew and xddisp depend on x because xdr depends on x
          dxdskew = (xdg+xddisp/(2.*xdr))*(1.+xdr**2)*(pixelwd/xdfl)
          hrnlin = -(abs(hrr)+1./(2.*abs(hrr)))*(pixelwd/hrfl)/2.
          spaci = nt
          xorder0 = xorder1-spacing/2.
        else
          xdskew = 2.*xdg*xdr+tan(slitrot+2.*krot)
          xdsmile = -xdr*pixelwd/xdfl
!	why isn''t xdsmile corrected here too?
          xdnlin = -(xdr+1./(2.*xdr))*(pixelwd/xdfl)/2.
        endif
        cosrot = cos(detrot)
        sinrot = sin(detrot)
        xmid = (nx+1)/2.
        ymid = (ny+1)/2.

        do iz = 1,nz
          do iy = 1,ny
            do ix = 1,nx
              b(ix,iy) = a(ix,iy,iz)
            enddo   
          enddo   
          do iy = 1,ny
            do ix = 1,nx
!	illum was based on a testtorted array, which isn't quite right here
!              if(illum(ix,iy)<=0) then
!                a(ix,iy,iz) = 0.
!                cycle
!              endif
              x = ix-xmid
              y = iy-ymid
              if(crossdisp) then
!	but illx should be right
                if(illx(ix)<=0) then
                  a(ix,iy,iz) = 0.
                  cycle
                endif
                order = (ix-0.5)/spaci+0.5
                iorder = nint(order)
!	distance from order center
                dorder = order-iorder
!	undistorted, but non-integral order spacing, coordinate
                x = xorder0+order*spacing-xmid
                dx = spacing*dorder
!	slit skewing within orders by echelon smile
                y = y+hrskew*dx+hrsmile*dx**2
!	non-linearity of echelon spectrum
!	(subtract ymid**2 so middle moves and ends stay put)
                y = y+hrnlin*(y**2-ymid**2)
!	skewing by cross dispersion, K mirror, and cross-dispersion smile
!	note: xd dispersion depends on linear y (wavelength)
                x = x+xddisp*(iy-ymid)+xdskew*y+dxdskew*x*y &
                    +xdsmile*y**2
!	non-linearity of xd spectrum
                x = x+xdnlin*(x**2-xmid**2)
              elseif(modeinst/=MCAMERA) then
                if(illx(iy)<=0) then
                  a(ix,iy,iz) = 0.
                  cycle
                endif
!	skewing by cross-dispersion smile
                x = x+xdskew*y+xdsmile*y**2
                x = x+xdnlin*(x**2-xmid**2)
              endif
!	barrel distortion
              barrel = 1.-brl*((x-x0brl)**2+(y-y0brl)**2)/xmid**3
              x = x*barrel
              y = y*barrel
!	array rotation
              u = x*cosrot-y*sinrot
              v = y*cosrot+x*sinrot
              if(sincwt) then
!	sinc interpolation
                iu0 = int(u+xmid)
                iv0 = int(v+ymid)
                iu1 = iu0+1
                iv1 = iv0+1
                if((iu0<3).or.(iu0>(nx-2)) &
                  .or.(iv0<3).or.(iv0>(ny-2)) &
                  .or.(b(iu0,iv0)*b(iu1,iv0)*b(iu0,iv1)*b(iu1,iv1) &
                    ==0.)) then
                  a(ix,iy,iz) = 0.
                  cycle
                endif
                ndv = min(iv0-1,(ny-iv1),8)
                idv = 0
                do jdv = 1,ndv
                  jv1 = iv1-jdv
                  jv2 = iv0+jdv
                  if((jv1<1).or.(jv2>ny) &
                    .or.(b(iu0,jv1)*b(iu0,jv2)==0.)) &
                      exit
                  idv = jdv
                enddo   
                ndu = min(iu0-1,(nx-iu1),8)
                idu = 0
                do jdu = 1,ndu
                  ju1 = iu1-jdu
                  ju2 = iu0+jdu
                  if((ju1<1).or.(ju2>nx) &
                    .or.(b(ju1,iv0)*b(ju2,iv0)==0.)) &
                      exit
                  idu = jdu
                enddo   
                if((idu<2).or.(idv<2)) then
                  a(ix,iy,iz) = 0.
                  cycle
                endif
                du0 = u+xmid-iu0
                dv0 = v+ymid-iv0
                sindu = sin(PI*du0)
                sindv = sin(PI*dv0)
                iu1 = max(iu1-idu,1)
                iu2 = min(iu0+idu,nx)
                iv1 = max(iv1-idv,1)
                iv2 = min(iv0+idv,ny)
                sumb = 0.
                sumw = 0.
                sindv = -(-1.)**(iv0-iv1)*abs(sindv)
                do iv = iv1,iv2
                  sindv = -sindv
                  dv = v+ymid-iv
                  if(abs(dv)<TINY) then
                    sincy = PI
                  else
                    sincy = sindv/dv
                  endif
                  if((iv==iv1).or.(iv==iv2)) sincy = sincy/2.
                  sindu = -(-1.)**(iu0-iu1)*abs(sindu)
                  do iu = iu1,iu2
                    sindu = -sindu
                    du = u+xmid-iu
                    if(abs(du)<TINY) then
                      sinc = sincy*PI
                    else
                      sinc = sincy*sindu/du
                    endif
                    if((iu==iu1).or.(iu==iu2)) sinc = sinc/2.
                    if(b(iu,iv)/=0.) then
                      sumb = sumb+sinc*b(iu,iv)
                      sumw = sumw+sinc
                    endif
                  enddo   
                enddo   
                if(sumw>0.) then
                  a(ix,iy,iz) = sumb/sumw
                else
                  a(ix,iy,iz) = 0.
                endif
              else
!	quadratic bell interpolation
                iu1 = int(u+xmid)-1
                iv1 = int(v+ymid)-1
                iu2 = iu1+1
                iv2 = iv1+1
                iu3 = iu1+2
                iv3 = iv1+2
                iu4 = iu1+3
                iv4 = iv1+3
                du = iu3-(u+xmid)
                dv = iv3-(v+ymid)
                if((iu1>1).and.(iu4<nx) &
                  .and.(iv1>1).and.(iv4<ny)) then
                  du2 = du**2
                  du12 = (1.-du)**2
                  dv2 = dv**2
                  dv12 = (1.-dv)**2
                  b11 = b(iu1,iv1)
                  b21 = b(iu1+1,iv1)
                  b31 = b(iu1+2,iv1)
                  b41 = b(iu1+3,iv1)
                  b12 = b(iu1,iv2)
                  b22 = b(iu1+1,iv2)
                  b32 = b(iu1+2,iv2)
                  b42 = b(iu1+3,iv2)
                  b13 = b(iu1,iv3)
                  b23 = b(iu1+1,iv3)
                  b33 = b(iu1+2,iv3)
                  b43 = b(iu1+3,iv3)
                  b14 = b(iu1,iv4)
                  b24 = b(iu1+1,iv4)
                  b34 = b(iu1+2,iv4)
                  b44 = b(iu1+3,iv4)
                  wu1 = 0.25*du2
                  wu2 = 0.5-0.25*du12
                  wu3 = 0.5-0.25*du2
                  wu4 = 0.25*du12
                  wv1 = 0.25*dv2
                  wv2 = 0.5-0.25*dv12
                  wv3 = 0.5-0.25*dv2
                  wv4 = 0.25*dv12
                  if((b22*b23*b32*b33)/=(0.)) then
                    a(ix,iy,iz) = &
                      wv1*(wu1*b11+wu2*b21+wu3*b31+wu4*b41) &
                      +wv2*(wu1*b12+wu2*b22+wu3*b32+wu4*b42) &
                      +wv3*(wu1*b13+wu2*b23+wu3*b33+wu4*b43) &
                      +wv4*(wu1*b14+wu2*b24+wu3*b34+wu4*b44)
                  else
                    a(ix,iy,iz) = 0.
                  endif
                elseif((iu2>1).and.(iu3<nx) &
                  .and.(iv2>1).and.(iv3<nx)) then
                  b22 = b(iu1+1,iv2)
                  b32 = b(iu1+2,iv2)
                  b23 = b(iu1+1,iv3)
                  b33 = b(iu1+2,iv3)
                  wu2 = du
                  wu3 = 1.-du
                  wv2 = dv
                  wv3 = 1.-dv
                  if((b22*b23*b32*b33)/=(0.)) then
                    a(ix,iy,iz) = &
                      wv2*(wu2*b22+wu3*b32)+wv3*(wu2*b23+wu3*b33)
                  else
                    a(ix,iy,iz) = 0.
                  endif
                else
                  a(ix,iy,iz) = 0.
                endif
              endif
            enddo   
          enddo   
          if(verbose.and.(modeobs<MSCAN)) then
            call grey(mx,nx,ny,a(1,1,iz),0.,0.,iwin1)
            call waitasec(pause,.true.,a(1,1,iz),mx,nx,ny)
          endif
        enddo   

        inquire(unit=iuredh,opened=hdopen)
        if(hdopen.and.(nz>1).and.(modeobs<MSCAN)) then
          write(iuredh,'("tort    = T")')
          comment = 'distortions have been corrected'
          if(dofits) call fithlog('TORT    ',.true.,comment,iufith)
        endif

        return
      end


      subroutine qtort(a,nz,ierr)
!	a quick version of tort for debounce

        use dims
        use modes
        use consts

        real a(mx,my,1),b(mx,my)
        real nodpa,lores,kmirror,krot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical good

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /iwin/ iwin1,iwin2

        if(modetort==MNONE) return

!	calculate mapping of undistorted (x,y) onto arr (u,v)
!	then interpolate from u,v onto x,y

        if(crossdisp) then
          hrslrot = -slitrot+krot
!	include echelon smile and slit rotation here
          hrskew = 2.*hrg*abs(hrr)+tan(hrslrot)
!	include xd smile and spectrum rotation by k mirror here
          xdskew = 2.*xdg*xdr+tan(krot)
          xdsmile = -xdr*pixelwd/xdfl
          xddisp = (xdr*xdfl)/(abs(hrr)*hrfl)
          xdnlin = -(xdr+1./(2.*xdr))*(pixelwd/xdfl)/2.
          hrnlin = -(abs(hrr)+1./(2.*abs(hrr)))*(pixelwd/hrfl)/2.
        else
          xdskew = 2.*xdg*xdr+tan(slitrot+2.*krot)
          xdsmile = -xdr*pixelwd/xdfl
          xdnlin = -(xdr+1./(2.*xdr))*(pixelwd/xdfl)/2.
        endif
        cosrot = cos(detrot)
        sinrot = sin(detrot)
        xmid = (nx+1)/2.
        ymid = (ny+1)/2.
        xorder0 = xorder1-spacing/2.

        do iz = 1,nz
          do iy = 1,ny
            do ix = 1,nx
              b(ix,iy) = a(ix,iy,iz)
            enddo   
          enddo   
          do iy = 1,ny
            do ix = 1,nx
              x = ix-xmid
              y = iy-ymid
              if(crossdisp) then
!	slit skewing within orders by echelon smile
                order = (ix-xorder0)/spacing
                iorder = nint(order)
                dorder = order-iorder
                dx = spacing*dorder
                y = y+hrskew*dx
!	non-linearity of echelon spectrum
                y = y+hrnlin*(y**2-ymid**2)
!	skewing by cross dispersion, k mirror, and cross-dispersion smile
!	note: xd dispersion depends on linear y (wavelength)
                x = x+xddisp*(iy-ymid)+xdskew*y+xdsmile*y**2
!	non-linearity of xd spectrum
                x = x+xdnlin*(x**2-xmid**2)
              elseif(modeinst/=MCAMERA) then
!	skewing by cross-dispersion smile
                x = x+xdskew*y+xdsmile*y**2
              endif
!	array rotation
              u = x*cosrot-y*sinrot
              v = y*cosrot+x*sinrot
!	just fetch the value from the nearest pixel
              iu = nint(u+xmid)
              iv = nint(v+ymid)
              if((iu>1).and.(iu<nx) &
                  .and.(iv>1).and.(iv<ny)) then
                a(ix,iy,iz) = b(iu,iv)
              else
                a(ix,iy,iz) = 0.
              endif
            enddo   
          enddo   
        enddo   

        return
      end


      subroutine testtort(c,ierr)
!	remove optical distortions from black-shiny to test tort,
!	find orders, and find illuminated pixels

        use ius
        use dims
        use modes
        use consts
        integer, parameter :: lx = 512
        integer, parameter :: ly = 512
        integer, parameter :: kx = 64
        integer, parameter :: ky = 64

        real b(mx,my),c(mx,my),bb(mx,my),db(mx,my),pow(mx),deriv(mx)
        real correl(4)
        real cr2(lx,ly),ci2(lx,ly),pow2(mx,my)
        real bx(mx),bbx(mx),dbx(mx),cx(mx),ccx(mx),bx1(mx),bx2(mx)
        real p(3),xm(3)
        real nodpa,lores,kmirror,krot
        character(60) object,feature,obsmode,instmode
        character(60) objtype,piname,pid,note,weather,warning,version
        logical good,baddata,doplot
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical pinhole,overlap,neworder
        character(8) yn

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrc/ object,feature,obsmode,instmode,objtype, &
                piname,pid,note,weather,warning,version
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
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /iwin/ iwin1,iwin2

        if((modetort==MNONE).or.(modetort==MOLD)) return
        doplot = (iwin1>0)
        if(modetort==MCELL) kfix = .true.
        nchange = 0

!	clean card array

        call setone(b,1,ierr)
        call clean(c,b,1,ierr)

!	tort card array (except echelon slit skewing)

        if(date<(1.00)) then
          pinhole = (slit<(180.))
        elseif(date<=(2.00)) then
          pinhole = (slit>(190.)).and.(slit<(290.))
        elseif(date<(2.08)) then
          pinhole = (slit>(250.)).and.(slit<(310.))
        else
          pinhole = (slit>(290.)).and.(slit<(350.))
        endif

  100   if(crossdisp) then
!	xdr could come from either echelle angle or wno0
!	wno0 isn''t really at array center but what can I do?
          xdskew = 2.*xdg*xdr+tan(krot)
          xdsmile = -xdr*pixelwd/xdfl
!	xd smile seems to be overcorrected.  put in fudge factor.
          xdsmile = 0.7*xdsmile
          xddisp = (xdr*xdfl)/(abs(hrr)*hrfl)
          xdnlin = -(xdr+1./(2.*xdr))*(pixelwd/xdfl)/2.
          dxdskew = (xdg+xddisp/(2.*xdr))*(1.+xdr**2)*(pixelwd/xdfl)
        elseif(modeinst/=MCAMERA) then
          xdskew = 2.*xdg*xdr+tan(slitrot+2.*krot)
          xdsmile = -xdr*pixelwd/xdfl
          xdnlin = -(xdr+1./(2.*xdr))*(pixelwd/xdfl)/2.
        endif
        cosrot = cos(detrot)
        sinrot = sin(detrot)
        xmid = (nx+1)/2.
        ymid = (ny+1)/2.
        do iy = 1,ny
          do ix = 1,nx
            x = ix-xmid
            y = iy-ymid
            if(crossdisp) then
              y = y+hrnlin*(y**2-ymid**2)
              x = x+xddisp*(iy-ymid)+xdskew*y+dxdskew*x*y+xdsmile*y**2
              x = x+xdnlin*(x**2-xmid**2)
            elseif(modeinst/=MCAMERA) then
              x = x+xdskew*y+xdsmile*y**2
            endif
            barrel = 1.-brl*((x-x0brl)**2+(y-y0brl)**2)/xmid**3
            x = x*barrel
            y = y*barrel
            u = x*cosrot-y*sinrot
            v = y*cosrot+x*sinrot
            iu1 = u+xmid
            iu2 = iu1+1
            iv1 = v+ymid
            iv2 = iv1+1
            du1 = u+xmid-iu1
            du2 = 1.0-du1
            dv1 = v+ymid-iv1
            dv2 = 1.0-dv1
            if((iu1>=1).and.(iu2<=nx) &
                .and.(iv1>=1).and.(iv2<=ny)) then
              b(ix,iy) = du2*dv2*c(iu1,iv1)+du1*dv2*c(iu2,iv1) &
                +du2*dv1*c(iu1,iv2)+du1*dv1*c(iu2,iv2)
!	this causes problems if used in tort
!	since testtorted array differs from torted
              illum(ix,iy) = 1
            else
              b(ix,iy) = 0.
              illum(ix,iy) = -1
            endif
          enddo   
        enddo   
        if(verbose.or.testrun) then
          print '(" Plotting testtorted array")'
          call grey(mx,nx,ny,b,0.,0.,iwin1)
          call waitasec(testrun,.true.,b,mx,nx,ny)
        endif

      if(crossdisp) then

!	db = db/dx
!	bb = (db/dx)**2

        ymid = (ny+1)/2.
        ywid = ny
        do iy = 1,ny
          tapery = 1.0-sin((iy-ymid)/ywid)**4
          kix = 0
          lix = nx
          do ix = 1,nx
            if(illum(ix,iy)>0) then
              if(kix==0) kix = ix
              lix = ix
            else
              bb(ix,iy) = 0.
              db(ix,iy) = 0.
            endif
          enddo
          xmid = (kix+lix)/2.
          xwid = (lix-kix)/PI
          do ix = kix,lix
            if((illum(ix+1,iy)>0).and.(illum(ix-1,iy)>0)) then
              taper = tapery*(1.0-sin((ix-xmid)/xwid)**4)
              bb(ix,iy) = taper*(b(ix+1,iy)-b(ix-1,iy))**2
              db(ix,iy) = taper*(b(ix+1,iy)-b(ix-1,iy))
            else
              bb(ix,iy) = 0.
              db(ix,iy) = 0.
            endif
          enddo   
        enddo   
        if(verbose) then
          print '(" Plotting black")'
          call grey(mx,nx,ny,b,0.,0.,iwin1)
          call waitasec(pause,.true.,b,mx,nx,ny)
          print '(" Plotting black derivative")'
          call grey(mx,nx,ny,db,0.,0.,iwin1)
          call waitasec(pause,.true.,db,mx,nx,ny)
        endif

        dw = 0.5/(sqrt(hrr**2/(1.+hrr**2))*hrdgr)
        predict = 2.0*xdr*xdfl*dw/(waveno0*pixelwd)
!	adjust predict to fit observed better
        if(date>1.0300) predict = 0.995*predict
        predi = 512./predict+1.

  400   if(ffttort) then
!	use 2-d fft to determine orientation of orders

          do ix = 1,lx
            do iy = 1,ly
              cr2(ix,iy) = 0.
              ci2(ix,iy) = 0.
            enddo   
          enddo   
          xnx = nx
          sumsq = 0.
          do iy = 1,ny
            do ix = 1,nx
              cr2(ix,iy) = bb(ix,iy)
              sumsq = sumsq+cr2(ix,iy)**2
            enddo   
          enddo   

          nux = 9
          nuy = 9
          call tdrfft(cr2,ci2,nux,nuy)

!	look for peak near 2nd harmonic of bb
!	would fundamental or db work better? 

          if(predi>20.) then
            ix1 = 1.8*predi
            ix2 = 2.2*predi+1.0
            iy2 = 0.2*predi+1.0
          else
            ix1 = 1.6*predi
            ix2 = 2.4*predi+1.0
            iy2 = 0.4*predi+1.0
          endif
          
          powmax = 0.
          powmax2 = 0.
          ixmax = 0
          iymax = 0
          xmax = 0.
          ymax = 0.
          angle = 0.
          ky2 = ky/2
          do iy = 1,ny
            do ix = 1,nx
              pow2(ix,iy) = 1.
            enddo
          enddo
          do iy = 1,iy2
            iyn = ly-iy+1
            do ix = ix1,ix2
              powi = cr2(ix,iy)**2+ci2(ix,iy)**2
              pow2(ix,iy+ky2) = powi
              if(powi>powmax) then
                ixmax = ix
                iymax = iy+ky2
                powmax = powi
              endif
              powi = cr2(ix,iyn)**2+ci2(ix,iyn)**2
              pow2(ix,ky2-iy+1) = powi
              if(powi>powmax) then
                ixmax = ix
                iymax = ky2-iy+1
                powmax = powi
              endif
            enddo   
          enddo   
          if(verbose) then
            print '(" Plotting 2D FFT power.  max =",es10.2,2i4)', &
              powmax,ixmax,iymax
            call grey(mx,nx/2,ny/2,pow2,0.,0.,iwin1)
            call waitasec(.true.,.true.,pow2,mx,nx/2,ny/2)
          endif
          if((ixmax==ix1).or.(ixmax==ix2) &
            .or.(iymax==(ky2-iy2+1)).or.(iymax==(ky2+iy2))) then
            print '(" Can''t find 2nd harmonic max in 2D FFT (1)")'
            print '(" ixmax,iymax = ",6i4)', &
              ixmax,iymax,ix1,ix2,ky2-iy2+1,ky2+iy2
            print '(8es10.2)', (pow2(ix,iymax),ix=ix1,ix2)
            goto 410
          endif
!	find peak in x for each y
          do jy = 1,3
            iy = iymax+jy-2
            p1 = pow2(ixmax-1,iy)
            p2 = pow2(ixmax,iy)
            p3 = pow2(ixmax+1,iy)
            pa = p2
            pb = (p3-p1)/2.
            pc = p1-2.*p2+p3
            if(pc>=0.) then
              print '(" Can''t find 2nd harmonic max in 2D FFT (2)")'
              print '(" ixmax,iy = ",2i4)', ixmax,iy
              print '(" p1,p2,p3 = ",3e10.3)', p1,p2,p3
              xmax = ixmax
              ymax = iymax
              goto 410
            endif
            xmax = -pb/pc
            xm(jy) = xmax
            p(jy) = pa+pb*xmax+(pc/2.)*xmax**2
          enddo   
!	find peak in y at peak in x
          pa = p(2)
          pb = (p(3)-p(1))/2.
          pc = p(1)-2.*p(2)+p(3)
          if(pc>=0.) then
            print '(" Can''t find 2nd harmonic max in 2D FFT (3)")'
              print '(" ixmax,iymax = ",2i4)', ixmax,iymax
              print '(" p1,p2,p3 = ",3e10.3)', p(1),p(2),p(3)
            xmax = ixmax
            ymax = iymax
            goto 410
          endif
          ymax = -pb/pc
          pmax2 = pa+pb*ymax+(pc/2.)*ymax**2
          powmax2 = powmax
          xa = xm(2)
          xb = (xm(3)-xm(1))/2.
          xc = xm(1)-2.*xm(2)+xm(3)
          xmax = xa+xb*ymax+xc*ymax**2+(ixmax-1)
          ymax = ymax-ky2+(iymax-1)
          spac = 512./(sqrt(xmax**2+ymax**2))
          if(verbose) &
            print '(" 2D FFT 2nd harmonic vector: ",2f10.6)', xmax,ymax
  410     xmax2 = xmax
          ymax2 = ymax
          angle = atan2(ymax,xmax)

          spacold = spacing
          noldt = int(spacold)

!	look for fundamental peak to find order spacing
          if(predi>64.) then
            ix1 = 0.9*predi
            ix2 = 1.1*predi+1.0
            iy2 = 0.1*predi+1.0
          elseif(predi>16.) then
            ix1 = 0.8*predi
            ix2 = 1.2*predi+1.0
            iy2 = 0.2*predi+1.0
          else
            ix1 = 0.6*predi
            ix2 = 1.4*predi+1.0
            iy2 = 0.4*predi+1.0
          endif
          powmax = 0.
          ixmax = 0
          iymax = 0
          ky2 = ky/2
          do iy = 1,iy2
            iyn = ly-iy+1
            do ix = ix1,ix2
              powi = cr2(ix,iy)**2+ci2(ix,iy)**2
              pow2(ix,iy+ky2) = powi
              if(powi>powmax) then
                ixmax = ix
                iymax = iy+ky2
                powmax = powi
              endif
              powi = cr2(ix,iyn)**2+ci2(ix,iyn)**2
              pow2(ix,ky2-iy+1) = powi
              if(powi>powmax) then
                ixmax = ix
                iymax = ky2-iy+1
                powmax = powi
              endif
            enddo   
          enddo   
          if((ixmax==ix1).or.(ixmax==ix2) &
            .or.(iymax==(ky2-iy2+1)).or.(iymax==(ky2+iy2))) then
            print '(" Can''t find fundamental max in 2D FFT (1)")'
            print '(" ixmax,iymax = ",2i4)', ixmax,iymax
            powmax = 0.
            goto 490
          endif
          do jy = 1,3
            iy = iymax+jy-2
            p1 = pow2(ixmax-1,iy)
            p2 = pow2(ixmax,iy)
            p3 = pow2(ixmax+1,iy)
            pa = p2
            pb = (p3-p1)/2.
            pc = p1-2.*p2+p3
            if(pc>=0.) then
              print '(" Can''t find fundamental max in 2D FFT (2)")'
              print '(" ixmax,iy = ",2i4)', ixmax,iy
              print '(" p1,p2,p3 = ",3e10.3)', p1,p2,p3
              powmax = 0.
              goto 490
            endif
            xmax1 = -pb/pc
            xm(jy) = xmax1
            p(jy) = pa+pb*xmax1+(pc/2.)*xmax1**2
          enddo   
          pa = p(2)
          pb = (p(3)-p(1))/2.
          pc = p(1)-2.*p(2)+p(3)
          if(pc>=0.) then
            print '(" Can''t find fundamental max in 2D FFT (3)")'
            print '(" ixmax,iymax = ",2i4)', ixmax,iymax
            print '(" p1,p2,p3 = ",3e10.3)', p(1),p(2),p(3)
            powmax = 0.
            goto 490
          endif
          ymax = -pb/pc
          pmax = pa+pb*ymax+(pc/2.)*ymax**2
          xa = xm(2)
          xb = (xm(3)-xm(1))/2.
          xc = xm(1)-2.*xm(2)+xm(3)
          xmax = xa+xb*ymax+(xc/2.)*ymax**2+(ixmax-1)
          ymax = ymax-ky2+(iymax-1)
          if(verbose) &
            print '(" 2D FFT fundamental vector: ",2f10.6)', xmax,ymax

  490     if(powmax2>powmax) then
            print '(" Using 2nd harmonic for order spacing")'
            xmax = xmax2/2.
            ymax = ymax2/2.
          elseif(powmax==0.) then
            print '("*** Can''t determine order spacing")'
            write(iupl,'("!!! Can''t determine order spacing")')
            print '(" Plotting testtorted array")'
            call grey(mx,nx,ny,b,0.,0.,iwin1)
            call waitasec(.true.,.true.,b,mx,nx,ny)
            print '(" Plotting 2D FFT of derivative")'
            print '(" predi,ix1,ix2,iy2,ky2 = ",f8.1,4i8)', &
              predi,ix1,ix2,iy2,ky2
            call grey(mx,ky,ky,pow2,0.,0.,iwin1)
            call waitasec(pause,.true.,pow2,mx,ky,ky)
            print '(" Hit RETURN to continue",$)'
            read '(a8)', yn
            ierr = 9
            return
          endif

          spac = lx/(sqrt(xmax**2+ymax**2))
          print '(" Order spacing, angle: ",2f8.3)', spac,angle
          print '(" Predicted order spacing, krot:",2f8.3)', predict,krot

          if(abs(spac-predict)>predict/50.) then
            print '("*** Order spacing disagrees with prediction")'
            if(abs(spac-predict)>predict/25.) then
              print '(" Hit RETURN to continue",$)'
              read '(a8)', yn
            endif
            write(iupl,'(" Order spacing, angle: ",2f8.3)') spac,angle
            write(iupl,'(" Predicted order spacing:",f8.3)') predict
            write(iupl,'(" order spacing disagrees with prediction")')
          endif

        else
!	get order orientation from correlation of top and bottom of orders

          ny2 = ny/2
          ny4 = ny/4
          do ix = 1,nx
            sumb = 0.
            sumn = 0.
            do iy = 9,ny2
              if(b(ix,iy)/=0.) then
                sumb = sumb+b(ix,iy)
                sumn = sumn+1.
              endif
              if(sumn>ny4) then
                bx1(ix) = sumb/sumn
              else
                bx1(ix) = 0.
              endif
            enddo
            sumb = 0.
            sumn = 0.
            do iy = ny2+1,ny-8
              if(b(ix,iy)/=0.) then
                sumb = sumb+b(ix,iy)
                sumn = sumn+1.
              endif
              if(sumn>ny4) then
                bx2(ix) = sumb/sumn
              else
                bx2(ix) = 0.
              endif
            enddo
          enddo
          ishmx = -5
          sumbmx = 0.
          do ishft = -4,4
            sumb = 0.
            do ix = 9,nx-8
              sumb = sumb+bx1(ix)*bx2(ix+ishft)
            enddo
            if(sumb>sumbmx) then
              ishmx = ishft
              sumbmx = sumb
            endif
          enddo
          do ishft = ishmx-1,ishmx+1
            sumb = 0.
            do ix = 9,nx-8
              sumb = sumb+bx1(ix)*bx2(ix+ishft)
            enddo
            correl(ishft-ishmx+2) = sumb
          enddo
          dcdsh = (correl(3)-correl(1))/2.
          ddcdsh = correl(1)-2.*correl(2)+correl(3)
          if(ddcdsh>abs(dcdsh)) then
            dsh = sign(1.,dcdsh)
          else
            dsh = -dcdsh/ddcdsh
          endif
          angle = -(ishmx+dsh)/120.
          if(abs(ishmx)>=4) &
            print '("*** Large order angle.",f8.3, &
            "  Set krot manually?",f8.3)', angle

!	get order spacing from correlation of b with shifted b or db/dx

          xny = ny
          bxsum = 0.
          ny8 = ny/8
          do ix = 1,nx
            bsum = 0.
            dbsum = 0.
            do iy = 1+ny8,ny-ny8
              bsum = bsum+b(ix,iy)
              dbsum = dbsum+db(ix,iy)
            enddo
            bx(ix) = bsum
            bbx(ix) = dbsum
            bxsum = bxsum+bsum
          enddo
          bxavg = bxsum/nx
          do ix = 1,nx
            bx(ix) = bx(ix)-bxavg
          enddo

          nx2 = nx/2
          nx8 = nx/8
          nx3 = 3*nx8
          xnx = nx
          if(predict<xnx/10.) then
            kdx = 3.5*predict
            ldx = 4.5*predict
            ldx = min(ldx,nx)
            bxmax = 0.
            bbxmax = 0.
            idxmax = 0
            iddxmax = 0
            do idx = kdx,ldx
              bxsum = 0.
              bbxsum = 0.
              do ix = nx8,nx3
                bxsum = bxsum+bx(ix)*bx(ix+idx)
                bbxsum = bbxsum+bbx(ix)*bbx(ix+idx)
              enddo
              cx(idx) = bxsum
              ccx(idx) = bbxsum
              if(bxsum>bxmax) then
                idxmax = idx
                bxmax = bxsum
              endif
              if(bbxsum>bbxmax) then
                iddxmax = idx
                bbxmax = bbxsum
              endif
            enddo
            if((idxmax>kdx).and.(idxmax<ldx)) then
              cx1 = cx(idxmax-1)
              cx2 = cx(idxmax)
              cx3 = cx(idxmax+1)
              dxmax = idxmax+0.5*(cx3-cx1)/(2.*cx2-(cx1+cx3))
            else
              dxmax = -idxmax
              print '("*** Error finding max b(x)",i6,es10.2)', idxmax,bxmax
              print '(" bx(ix)")'
              print '(8es10.2)', (bx(ix),ix=1,nx)
              print '(" cx(idx)")'
              print '(8es10.2)', (cx(idx),idx=kdx,ldx)
            endif
            if((iddxmax>kdx).and.(iddxmax<ldx)) then
              ccx1 = ccx(iddxmax-1)
              ccx2 = ccx(iddxmax)
              ccx3 = ccx(iddxmax+1)
              ddxmax = iddxmax+0.5*(ccx3-ccx1)/(2.*ccx2-(ccx1+ccx3))
            else
              ddxmax = -iddxmax
              print '(" Can''t find correlation peak in db/dx")'
            endif
            dxmax = dxmax/4.
            ddxmax = ddxmax/4.
          elseif(predict<xnx/5.) then
            kdx = 1.5*predict
            ldx = 2.5*predict
            ldx = min(ldx,nx)
            bxmax = 0.
            bbxmax = 0.
            idxmax = 0
            iddxmax = 0
            do idx = kdx,ldx
              bxsum = 0.
              bbxsum = 0.
              do ix = nx8,nx3
                bxsum = bxsum+bx(ix)*bx(ix+idx)
                bbxsum = bbxsum+bbx(ix)*bbx(ix+idx)
              enddo
              cx(idx) = bxsum
              ccx(idx) = bbxsum
              if(bxsum>bxmax) then
                idxmax = idx
                bxmax = bxsum
              endif
              if(bbxsum>bbxmax) then
                iddxmax = idx
                bbxmax = bbxsum
              endif
            enddo
            if((idxmax>kdx).and.(idxmax<ldx)) then
              cx1 = cx(idxmax-1)
              cx2 = cx(idxmax)
              cx3 = cx(idxmax+1)
              dxmax = idxmax+0.5*(cx3-cx1)/(2.*cx2-(cx1+cx3))
            else
              dxmax = -idxmax
              print '("*** Error finding max b(x)",i6,es10.2)', idxmax,bxmax
              print '(" bx(ix)")'
              print '(8es10.2)', (bx(ix),ix=1,nx)
              print '(" cx(idx)")'
              print '(8es10.2)', (cx(idx),idx=mdx,ndx)
            endif
            if((iddxmax>kdx).and.(iddxmax<ldx)) then
              ccx1 = ccx(iddxmax-1)
              ccx2 = ccx(iddxmax)
              ccx3 = ccx(iddxmax+1)
              ddxmax = iddxmax+0.5*(ccx3-ccx1)/(2.*ccx2-(ccx1+ccx3))
            else
              ddxmax = -iddxmax
              print '(" Can''t find correlation peak in db/dx")'
            endif
            dxmax = dxmax/2.
            ddxmax = ddxmax/2.
          else
            kdx = 0.5*predict
            ldx = 1.5*predict
            ldx = min(ldx,nx)
            bxmax = 0.
            bbxmax = 0.
            idxmax = 0
            iddxmax = 0
            do idx = kdx,ldx
              bxsum = 0.
              bbxsum = 0.
              do ix = nx8/2,nx2
                bxsum = bxsum+bx(ix)*bx(ix+idx)
                bbxsum = bbxsum+bbx(ix)*bbx(ix+idx)
              enddo
              cx(idx) = bxsum
              ccx(idx) = bbxsum
              if(bxsum>bxmax) then
                idxmax = idx
                bxmax = bxsum
              endif
              if(bbxsum>bbxmax) then
                iddxmax = idx
                bbxmax = bbxsum
              endif
            enddo
            if((idxmax>kdx).and.(idxmax<ldx)) then
              cx1 = cx(idxmax-1)
              cx2 = cx(idxmax)
              cx3 = cx(idxmax+1)
              dxmax = idxmax+0.5*(cx3-cx1)/(2.*cx2-(cx1+cx3))
            else
              dxmax = -idxmax
              print '("*** Error finding max b(x)",i6,es10.2)', idxmax,bxmax
              print '(" bx(ix)")'
              print '(8es10.2)', (bx(ix),ix=1,nx)
              print '(" cx(idx)")'
              print '(8es10.2)', (cx(idx),idx=mdx,ndx)
            endif
            if((iddxmax>kdx).and.(iddxmax<ldx)) then
              ccx1 = ccx(iddxmax-1)
              ccx2 = ccx(iddxmax)
              ccx3 = ccx(iddxmax+1)
              ddxmax = iddxmax+0.5*(ccx3-ccx1)/(2.*ccx2-(ccx1+ccx3))
            else
              ddxmax = -iddxmax
              print '(" Can''t find correlation peak in db/dx")'
            endif
          endif
          print '(" Order spacing from b and db/dx: ",2f8.3)', dxmax,ddxmax

          if((ddxmax>0.).and.(abs(ddxmax-predict)<predict/25.)) then
            spac = ddxmax
          elseif((dxmax>0.).and.(abs(dxmax-predict)<predict/25.)) then
            spac = dxmax
          else
            print '(" Predicted order spacing:",f8.3)', predict
            print '(" Enter order spacing or 0 to quit: ",$)'
            read *, spac
            if(spac==0.) stop
          endif
          print '(" Order spacing, angle: ",2f8.3)', spac,angle
          print '(" Predicted order spacing, krot:",2f8.3)', predict,krot
          if(abs(spac-predict)>predict/50.) then
            print '("*** Order spacing disagrees with prediction")'
            if(abs(spac-predict)>predict/25.) then
              print '(" Hit RETURN to continue",$)'
              read '(a8)', yn
            endif
            write(iupl,'(" Order spacing, angle: ",2f8.3)') spac,angle
            write(iupl,'(" Predicted order spacing:",f8.3)') predict
            write(iupl,'(" order spacing disagrees with prediction")')
          endif
        endif

        spacold = spacing
        spacing = spac
        noldt = int(spacold)
        nt = int(spacing)
        if(dosum.and.(nt/=noldt).and.(wno0>0.).and.(spacold/=32.)) then
          print '("*** Order spacing changed from",f8.3," to",f8.3)', &
            spacold,spacing
          if((abs(spacing/spacold-1.)<(0.001))) then
            print '(" Setting it back to avoid trashing sum")'
            spacing = spacold
            nt = noldt
          elseif((abs(spacing/spacold-1.)<(0.1))) then
            print '("*** Sum may be trashed")'
          endif
        endif

        if(kfix) then
          print '(" Leaving krot fixed")'
        elseif((abs(angle)>0.05).or.(doplot.and.(abs(angle)>0.02))) then
          print '("*** order angle > 0.02 ",f8.3)', angle
          write(iupl,'(" order angle > 0.02")')
          if(doplot) then
            print '(" Plotting testtorted array")'
            call grey(mx,nx,ny,b,0.,0.,iwin1)
          endif
          call waitasec(.true.,.true.,b,mx,nx,ny)
          print '(" Should I change krot to",f7.3"? (y) ",$)', krot-angle
          read '(a8)', yn
          if(index(yn,'n')==0) then
            krot = krot-angle
            print '(" Changing krot to ",f7.3)', krot
            write(iupl,'(" Changing krot to ",f7.3)') krot
            goto 100
          elseif(.not.ffttort) then
            print '(" Switch to FFT order finding? (n) ",$)'
            read '(a8)', yn
            if(index(yn,'y')>0) then
              ffttort = .true.
              goto 400
            endif
          endif
          print '(" fife may be confused.  Enter krot ",$)'
          read *, krot
          kfix = .true.
          goto 100
        elseif(abs(angle)>0.001) then
          krot = krot-angle
          print '(" Changing krot to correct angle; krot =",f7.3)', krot
          write(iupl,'(" Changing krot to krot =",f7.3)') krot
          nchange = nchange+1
          if(nchange<4) then
            goto 100
          elseif(.not.ffttort) then
            print '(" fife may be confused." &
              "  Switch to FFT order finding? (n) "$)'
            read '(a8)', yn
            if(index(yn,'y')>0) then
              ffttort = .true.
              goto 400
            endif
          endif
          print '(" fife is confused.  Enter krot ",$)'
          read *, krot
          kfix = .true.
          goto 100
        endif

!	new method to find xorder1
        neworder = .true.
!        neworder = .false.

        idbmx = 0
        idbmn = 0
        dbmx = 0.
        dbmn = 0.
        do ix = 1,nx
          bbsum = 0.
          dbsum = 0.
          nbsum = 0
          do iy = ny/4+1,3*ny/4
            if(bb(ix,iy)/=0.) then
              bbsum = bbsum+bb(ix,iy)
              dbsum = dbsum+db(ix,iy)
              nbsum = nbsum+1
            endif
          enddo
          nbsum = max(nbsum,1)
          bbxi = bbsum/nbsum
          dbxi = dbsum/nbsum
          bbx(ix) = bbxi
          dbx(ix) = dbxi
          if(dbxi>dbmx) then
            idbmx = ix
            dbmx = dbxi
          endif
          if(dbxi<dbmn) then
            idbmn = ix
            dbmn = dbxi
          endif
        enddo
        dbslope = (dbx(idbmn+1)-dbx(idbmn-1))/2.
        dbcurve = dbx(idbmn-1)+dbx(idbmn+1)-2.*dbx(idbmn)
        if(dbcurve>abs(dbslope)) then
          xdbmn = idbmn-dbslope/dbcurve
        else
          xdbmn = idbmn-sign(1.,dbslope)
        endif
        dbslope = (dbx(idbmx+1)-dbx(idbmx-1))/2.
        dbcurve = dbx(idbmx-1)+dbx(idbmx+1)-2.*dbx(idbmx)
        if(-dbcurve>abs(dbslope)) then
          xdbmx = idbmx-dbslope/dbcurve
        else
          xdbmx = idbmx+sign(1.,dbslope)
        endif
        xnx2 = nx/2.
        if(abs(xdbmn-xnx2)>spacing/2.) &
          xdbmn = xdbmn-spacing*nint((xdbmn-xnx2)/spacing)
        if(abs(xdbmx-xdbmn)>spacing/2.) &
          xdbmx = xdbmx-spacing*nint((xdbmx-xdbmn)/spacing)
        if((xdbmn>xdbmx).and.(.not.pinhole).and.(xdbmn-xdbmx<spacing/4.)) then
          print '(" Overlapping echelon orders")'
          overlap = .true.
          xorder0 = xdbmn-1.
        else
          overlap = .false.
          xorder0 = xdbmx-1.
        endif
        xorder0 = mod(xorder0,spacing)
        if(neworder) then
          if(xorder0>spacing-1.) xorder0 = mod(xorder0+1.,spacing)
          norder0 = (nx-xorder0+abs(xdbmx-xdbmn))/spacing
          xolder1 = xorder1
          nolder = norder
          xorder1 = xorder0
          norder = norder0
          ispac = spacing
          nxo = spacing+1.
          goto 310
        else
          norder0 = (nx-xorder0+abs(xdbmx-xdbmn)-2.)/spacing
        endif

!	sum derivative over y and orders to find orders

        nolder = norder
        xolder1 = xorder1
        nxo = spacing+1.
        norder = mx/nxo-1
        if(norder>2) then
          nxo2 = 2.*spacing+1.
          nx1 = nxo/2
          nx2 = nx1+nxo
        elseif(norder==2) then
          nxo2 = 2.*spacing+1.
          nx1 = 2
          nx2 = nxo
        elseif(norder==1) then
          nxo2 = spacing+1.
          nx1 = 2
          nx2 = nxo
        else
          print '(" I need a routine to find orders in such a short spectrum")'
          stop
        endif
        nmin = ny*norder/2
        do ixo = 1,nxo2
          sumb = 0.
          sumd = 0.
          nsumb = 0
          nsumd = 0
          do iorder = 1,norder
            ix = ixo+(iorder-1)*spacing
            ix1 = ix-1
            ix2 = ix+1
            if((ix1<1).or.(ix2>nx)) cycle
            do iy = 1,ny
              if(b(ix,iy)/=0.) then
                nsumb = nsumb+1
                sumb = sumb+b(ix,iy)
              endif
              if(b(ix1,iy)*b(ix2,iy)/=0.) then
                nsumd = nsumd+1
                sumd = sumd+b(ix2,iy)-b(ix1,iy)
              endif
            enddo   
          enddo   
          if(nsumb>nmin) then
            pow(ixo) = sumb/nsumb
          else
            pow(ixo) = 0.
          endif
          if(nsumd>nmin) then
            deriv(ixo) = sumd/nsumd
          else
            deriv(ixo) = 0.
          endif
        enddo   

        powmin = pow(nx1-1)
        powmax = pow(nx1-1)
        dermin = deriv(nx1-1)
        dermax = deriv(nx1-1)
        ifal = nx1-1
        iris = nx1-1
        do ixo = nx1,nx2
          if(pow(ixo)==0.) then
            cycle
          elseif(pow(ixo)<powmin) then
            powmin = pow(ixo)
          elseif(pow(ixo)>powmax) then
            powmax = pow(ixo)
          endif
          if(deriv(ixo)==0.) then
            cycle
          elseif(deriv(ixo)<dermin) then
            dermin = deriv(ixo)
            ifal = ixo
          elseif(deriv(ixo)>dermax) then
            dermax = deriv(ixo)
            iris = ixo
          endif
        enddo   
        der1 = deriv(ifal-1)
        der2 = deriv(ifal)
        der3 = deriv(ifal+1)
        xfal = ifal+0.5*(der1-der3)/(der1-2.*der2+der3)
        der1 = deriv(iris-1)
        der2 = deriv(iris)
        der3 = deriv(iris+1)
        xris = iris+0.5*(der1-der3)/(der1-2.*der2+der3)
        if(xfal>xris) xfal = xfal-spacing
        xorder1 = (xfal+xris)/2.

        if((.not.pinhole).and.((powmax/powmin)<(2.)) &
          .and.((xris-xfal)>(spacing/2.))) then
          print '(" ** Overlapping echelon orders",2f8.2)', xris,xfal
          write(iupl,'(" ** Overlapping echelon orders")')
          overlap = .true.
          xorder1 = xorder1-spacing/2.
        else
          overlap = .false.
        endif
!	for now make sure all orders are fully on array
        if(xorder1<0.) xorder1 = xorder1+spacing
        if(xorder1>spacing) xorder1 = xorder1-spacing
        norder = (nx-xorder1)/spacing
        if(norder*ny>ms) then
          print '("*** norder too big for spec array",2i4)', norder,ms/ny
          norder = ms/ny
        endif
        if(verbose) &
          print '(" xfal,xris,xorder1:",3f8.3)', xfal,xris,xorder1

  310   do ixo = 1,nxo
          sumb = 0.
          nsumb = 0
          do iorder = 1,norder
            ix = xorder1+(iorder-1)*spacing+ixo
            if((ix<0).or.(ix>nx)) cycle
            do iy = 1,ny
              if(b(ix,iy)/=0.) then
                sumb = sumb+b(ix,iy)
                nsumb = nsumb+1
              endif
            enddo   
          enddo   
          if(nsumb>nmin) then
            pow(ixo) = sumb/nsumb
          else
            pow(ixo) = 0.
          endif
        enddo   

        powmin = pow(1)
        powmax = pow(1)
        imin = 1
        imax = 1
        do ixo = 2,nxo
          if(pow(ixo)==0.) then
            cycle
          elseif((pow(ixo)<powmin).or.(powmin==0.)) then
            powmin = pow(ixo)
            imin = ixo
          elseif(pow(ixo)>powmax) then
            powmax = pow(ixo)
            imax = ixo
          endif
        enddo   
        if(overlap) then
          thresh = powmax-thrfac*(powmax-powmin)
        else
          thresh = powmin+thrfac*(powmax-powmin)
        endif
        if(neworder) goto 320

!	now let unilluminated edges of edge orders fall off

        ixmin = 1
        ixmax = nx
!	uncomment these lines to allow a fraction of each order to be lost
!       if((intext(1)>1).and.(intext(2)>intext(1)) &
!         .and.(intext(2)<nxo)) then
!         ixmin = 2-intext(1)
!         if((intext(3)>1).and.(intext(3)<intext(1))) &
!           ixmin = 2-intext(3)
!         ixmax = nx+nxo-intext(2)
!         if((intext(4)>intext(2)).and.(intext(4)<nxo)) &
!           ixmax = nx+nxo-intext(4)
!       endif
        if(xorder1>spacing/2.) xorder1 = xorder1-spacing
        do ixo = 1,nxo
          if((overlap.and.(pow(ixo)>thresh)) &
            .or.((.not.overlap).and.(pow(ixo)<thresh))) cycle
!	found first pixel above thresh; is it on the array?
          if(int(xorder1+ixo)<ixmin) xorder1 = xorder1+spacing
          if(int(xorder1+ixo)<1) &
            print '(" Kept first order although partially off array")'
          exit
        enddo   
        norder = nint((nx-xorder1)/spacing)
        ispac = int(spacing)
        do ixo = nxo,1,-1
          ix = (norder-1)*spacing+ixo
          if((overlap.and.(pow(ixo)>thresh)) &
            .or.((.not.overlap).and.(pow(ixo)<thresh))) cycle
!	found last pixel above thresh; is it on the array?
          if(ix>ixmax) then
            print '(" Pixel ",i2," of order ",i2," off array ",2i4)', &
              ixo,norder,ix,ixmax
            norder = norder-1
            print '(" Changing norder to ",i2)', norder
          endif
          exit
        enddo   

!	set illx = 0 where below threshold in x
!	set illx = -1 where outside of mostly covered orders

  320   ispac = int(spacing)
        spaci = ispac
        iorder1 = nint(xorder1-0.5)
        nbelow = 0
        noff = 0
        illsum = 0
        nox = (norder-1)*ispac
        if((intill(1)>0).and.(intill(2)<ispac) &
          .and.(intill(2)>intill(1))) then
          illx1 = intill(1)
          illx2 = intill(2)
        else
          illx1 = 0
          illx2 = ispac
        endif
        illxj = 0
        do ixo = 1,ispac
          if((ixo<illx1).or.(ixo>illx2)) then
            nbelow = nbelow+1
            illxi = 0
          elseif(overlap) then
            if(pow(ixo)>thresh) then
              nbelow = nbelow+1
              illxi = 0
            else
              illxi = 1
            endif
          else
            if(pow(ixo)<thresh) then
              nbelow = nbelow+1
              illxi = 0
            else
              illxi = 1
            endif
          endif
          if(ixo<(nxo/2)) then
            if((illxi==0).and.(illxj==1)) &
              print '("*** Weird illumination ",3i3)', ixo,illxi,illxj
          else
            if((illxi==1).and.(illxj==0)) &
              print '("*** Weird illumination ",3i3)', ixo,illxi,illxj
          endif
          do iorder = 1,norder
            ix = (iorder-1)*ispac+ixo
            if((ix>0).and.(ix<=nx)) then
              illx(ix) = illxi
              if(illxi==1) illsum = illsum+1
            elseif(illxi==1) then
              print '("*** Huh?",8i6)', iorder1,iorder,ispac,ixo,ix,illxi
            endif
          enddo   
          illxj = illxi
        enddo   
        ixn = norder*ispac+1
        do ix = ixn,nx
          illx(ix) = -1
        enddo   

        print '(" norder, spacing, xorder1, nt =", &
          i4,2f7.2,i4)', norder,spacing,xorder1,ispac
        write(iupl,'(" norder, spacing, xorder1, nt =", &
          i4,2f7.2,i4)') norder,spacing,xorder1,ispac
        if((norder/=nolder).and.(nolder>1).and.(wno0>0.)) then
          if(((abs(xorder1-xolder1)<(spacing/2.)) &
            .and.(mod((norder+nolder),4)==1)).or. &
            ((abs(xorder1-xolder1)>(spacing/2.)) &
            .and.(mod((norder+nolder),4)==3))) then
              print '("*** Warning: norder has changed")'
              write(iupl,'("!!! Warning: norder has changed")')
              print '(" nolder, xolder1 :",i3,f6.2)', nolder,xolder1
              if(dosum) print '("*** Leaving out of sum")'
              if(dosum) write(iupl,'("!!! Leaving out of sum")')
              wno0 = -abs(wno0)
              baddata = .true.
          else
            print '(" ** norder has changed but wno0 should be OK")'
            write(iupl, &
              '("!!! norder has changed but wno0 should be OK")')
          endif
        endif
        print '(" Power below threshold in ",i2," pixels per order")', &
          nbelow
        write(iupl,'(" Power below threshold in ",i2, &
          " pixels per order")') nbelow
        if((nbelow>(nxo/2)).and.(.not.pinhole)) then
          write(iupl,'("!!! That seems excessive.")')
          print '("*** That seems excessive.  pow(ixo):")'
          print '(8f9.2)', (pow(ixo),ixo=1,nxo),thresh
          print '(" Hit RETURN to continue",$)'
          read '(a8)', yn
        endif
        slitlen = spacing-nbelow
        if((bounce/=0.).and.(modeobs==MNOD)) then
          if((dist>(0.6*slitlen)).and.(dist<slitlen)) &
            print '("*** Warning: nodded source may confuse debounce")'
        endif
        nt = int(spacing)

      elseif(modeinst==MCAMERA) then

!	just collapse in y for camera mode

        powmin = 1e30
        powmax = -1e30
        sumb = 0.
        do ix = 1,nx
          do iy = 1,ny
            sumb = sumb + b(ix,iy)
          enddo   
          pow(ix) = sumb
          if(sumb<powmin) powmin = sumb
          if(sumb>powmax) powmax = sumb
        enddo   
        thresh = powmin+thrfac*(powmax-powmin)

        do ix = 1,nx
          if(pow(ix)<thresh) then
            illx(ix) = 0
          else
            illx(ix) = 1
          endif
        enddo   

      else

!	start with all illx = 1 for long-slit modes

        do ix = 1,nx
          illx(ix) = 1
        enddo   
        norder = 1
        ns = 256
        nt = 256

      endif

  500   illsum = 0
        do ix = 1,nx
          if(illx(ix)<1) then
            illxi = illx(ix)
            do iy = 1,ny
              illum(ix,iy) = min(illum(ix,iy),illxi)
            enddo   
          else
            illsum = illsum+1
          endif
        enddo   

!	determine illumination in y

        if(modeinst==MCAMERA) then
          illmin = nx/4
        else
          illmin = 3*illsum/4
        endif
        do iy = 1,ny
          sumb = 0.
          illsum = 0
          do ix = 1,nx
            if(illum(ix,iy)==1) then
              sumb = sumb+b(ix,iy)
              illsum = illsum+1
            endif
          enddo   
          if(illsum>illmin) then
            pow(iy) = sumb/illsum
          else
            pow(iy) = 0.
          endif
        enddo   

        powmax = pow(1)
        do iy = 1,ny
          if(pow(iy)>powmax) powmax = pow(iy)
        enddo   
        thresh = thrfac*powmax

!	illum may be wrong since it is based on testtorted array
        do iy = 1,ny
          if(pow(iy)<thresh) then
            do ix = 1,nx
              illum(ix,iy) = 0
            enddo
            if((modeinst==MLOW).or.(modeinst==MMED)) illx(iy) = 0
          else
            if((modeinst==MLOW).or.(modeinst==MMED)) illx(iy) = 1
          endif
        enddo

        return
      end


      subroutine submed(arr,flat,nz,ierr)
!	make slit median or quartile zero at each point in spectrum

        use ius
        use dims
        use modes

        real arr(mx,my,1),flat(mx,my),arriz(my,mx),arrmed(my)
        real krot
        logical good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical modeext

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /iwin/ iwin1,iwin2

        print '(" Subtracting median from each spectral point")'
        write(iupl,'(" Subtracting median from each spectral point")')

      if((modeinst==MHIMED).or.(modeinst==MHILOW)) then

        do iz = 1,nz
          do iorder = 1,norder
            ix0 = (iorder-1)*nt
            do iy = 1,ny
              do it = 1,nt
                ix = ix0+it
                arriz(iy,it) = arr(ix,iy,iz)
              enddo   
            enddo   
            if(nodon) then
              iq = 2
            else
              iq = 1
            endif
            call median(arrmed,arriz,1,my,nt,iq,ierr)
            do iy = 1,ny
              do it = 1,nt
                ix = ix0+it
                if(arr(ix,iy,iz)/=0.) &
                  arr(ix,iy,iz) = arr(ix,iy,iz)-arrmed(iy)
              enddo   
            enddo   
          enddo   
        enddo   

      elseif((modeinst==MMED).or.(modeinst==MLOW)) then

        if((intsky(1)>0).and.(intsky(2)>intsky(1))) then
          print '(" Taking median from skyint pixels ",2i4)', iy1,iy2
          iy1 = intsky(1)
          iy2 = min(ny,intsky(2))
        else
          iy1 = 1
          iy2 = ny
        endif
        iy0 = iy1-1
        ny0 = iy2-iy0
        do iz = 1,nz
          do ix = 1,nx
            do iy = 1,ny0
              arriz(ix,iy) = arr(ix,iy+iy0,iz)
            enddo   
            if(nodon) then
              iq = 2
            else
              iq = 1
            endif
            call median(arrmed,arriz,1,mx,ny0,iq,ierr)
            do iy = 1,ny
              if(illum(ix,iy)==1) then
                arr(ix,iy,iz) = arr(ix,iy,iz)-arrmed(ix)
              endif
            enddo   
          enddo   
        enddo   

      else

        do iz = 1,nz
          sum = 0.
          flatsum = 0.
          do iy = 1,ny
            do ix = 1,nx
              flt = flat(ix,iy)
              if((illum(ix,iy)==1).and.(flt>0.)) then
                sum = sum+arr(ix,iy,iz)/flt**2
                flatsum = flatsum+1./flt**2
              endif
            enddo   
          enddo   
          if(flatsum==0.) cycle
          avg = sum/flatsum
          do iy = 1,ny
            do ix = 1,nx
              if(illum(ix,iy)==1) then
                arr(ix,iy,iz) = arr(ix,iy,iz)-avg
              endif
            enddo   
          enddo   
        enddo   

      endif

        ierr = 0
        return
      end


      subroutine submean(arr,flat,nz,ierr)
!	make slit mean zero at each point in spectrum to remove sky noise
!	can call with std or flat for weighting?

        use ius
        use dims
        use modes

        real arr(mx,my,1),flat(mx,my)
        real krot
        logical good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical modeext

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /iwin/ iwin1,iwin2

        print '(" Subtracting mean from each spectral point")'
        write(iupl,'(" Subtracting mean from each spectral point")')

        if((modeinst==MHIMED).or.(modeinst==MHILOW)) then

        do iz = 1,nz
          do iorder = 1,norder
            ix0 = (iorder-1)*nt
            do iy = 1,ny
              sum = 0.
              flatsum = 0.
              do it = 1,nt
                ix = ix0+it
                if((ix>=1).and.(ix<=nx)) then
                  flt = flat(ix,iy)
                  if((illum(ix,iy)==1).and.(flt>0.)) then
                    sum = sum+arr(ix,iy,iz)/flt**2
                    flatsum = flatsum+1./flt**2
                  endif
                endif
              enddo   
              if(flatsum==0.) cycle
              avg = sum/flatsum
              do it = 1,nt
                ix = ix0+it
                if((ix>=1).and.(ix<=nx)) then
                  if(illum(ix,iy)==1) then
                    arr(ix,iy,iz) = arr(ix,iy,iz)-avg
                  endif
                endif
              enddo   
            enddo   
          enddo   
        enddo   

        elseif((modeinst==MMED).or.(modeinst==MLOW)) then

        if((intsky(1)>0).and.(intsky(2)>intsky(1))) then
          iy1 = intsky(1)
          iy2 = min(ny,intsky(2))
        else
          iy1 = 1
          iy2 = ny
        endif
        do iz = 1,nz
          do ix = 1,nx
            sum = 0.
            flatsum = 0.
            do iy = iy1,iy2
              flt = flat(ix,iy)
              if((illum(ix,iy)==1).and.(flt>0.)) then
                sum = sum+arr(ix,iy,iz)/flt**2
                flatsum = flatsum+1./flt**2
              endif
            enddo   
            if(flatsum==0.) cycle
            avg = sum/flatsum
            do iy = 1,ny
              if(illum(ix,iy)==1) then
                arr(ix,iy,iz) = arr(ix,iy,iz)-avg
              endif
            enddo   
          enddo   
        enddo   

        else

        do iz = 1,nz
          sum = 0.
          flatsum = 0.
          do iy = 1,ny
            do ix = 1,nx
              flt = flat(ix,iy)
              if((illum(ix,iy)==1).and.(flt>0.)) then
                sum = sum+arr(ix,iy,iz)/flt**2
                flatsum = flatsum+1./flt**2
              endif
            enddo   
          enddo   
          if(flatsum==0.) cycle
          avg = sum/flatsum
          do iy = 1,ny
            do ix = 1,nx
              if(illum(ix,iy)==1) then
                arr(ix,iy,iz) = arr(ix,iy,iz)-avg
              endif
            enddo   
          enddo   
        enddo   

        endif

        ierr = 0
        return
      end


      subroutine subcorr(arr,brr,flat,nz,ierr)
!	remove sky noise by correlating each pixel of brr with frame mean
!	or mean over skyint
!	and subtracting correlation array to make correlation zero
!	used for map and nod-off mode data

        use ius
        use dims
        use modes

        real arr(mx,my,1),brr(mx,my,1),flat(mx,my)
        real bavg(mx,my),dbdz(mx,my),corr(mx,my),avga(mp),avgb(mp)
        real plot(mx,my)
        real krot
        logical good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical scan,corrall

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /iwin/ iwin1,iwin2

        if(nz<4) then
          print '(" Can''t subtract sky from less than four frames")'
          return
        endif
        print '(" Subtracting correlated sky fluctuations in subcorr")'
        print '("*** This routine may subtract continuum")'
        print '(" Set skyint to specify raw pixels to correlate with")'
        write(iupl,'(" Subtracting sky fluctuations in subcorr")')
        scan = (modeobs==MSCAN).or.(modeobs==MMAP)

        ix1 = intsky(1)
        ix2 = intsky(2)
        iy1 = intsky(3)
        iy2 = intsky(4)
        corrall = ((ix1<1).or.(ix1>ix2).or.(ix2>nx) &
          .or.(iy1<1).or.(iy1>iy2).or.(iy2>ny))
        if(corrall) then
          nxy = nx*ny
        else
          nxy = (ix2-ix1+1)*(iy2-iy1+1)
        endif

!	calculate means and trends

        z2 = (nz+1)/2.
        if(.not.scan) then
          ng = 0
          z2 = 0.
          do iz = 1,nz
            if(wt(iz)/=0.) then
              z2 = z2+iz
              ng = ng+1
            endif
          enddo   
          if(ng>0) then
            z2 = z2/ng
          else
            print '("*** All weights are zero.  Abort.")'
            return
          endif
        endif
        z12 = 0.
        zavg = 0.
        davgdz = 0.
        ng = 0
        do iy = 1,ny
          do ix = 1,nx
            bavg(ix,iy) = 0.
            dbdz(ix,iy) = 0.
            corr(ix,iy) = 0.
          enddo   
        enddo   
        do iz = 1,nz
          if((wt(iz)==0.).and.(.not.scan)) cycle
          ng = ng+1
          suma = 0.
          sumb = 0.
          do iy = 1,ny
            do ix = 1,nx
              bavg(ix,iy) = bavg(ix,iy)+brr(ix,iy,iz)
              dbdz(ix,iy) = dbdz(ix,iy)+(iz-z2)*brr(ix,iy,iz)
              if(corrall.or.((ix>=ix1).and.(ix<=ix2) &
                .and.(iy>=iy1).and.(iy<=iy2))) then
                suma = suma+arr(ix,iy,iz)
                sumb = sumb+brr(ix,iy,iz)
              endif
            enddo   
          enddo   
          avga(iz) = suma/nxy
          avgb(iz) = sumb/nxy
          zavg = zavg+avgb(iz)
          davgdz = davgdz+(iz-z2)*avgb(iz)
          z12 = z12+(iz-z2)**2
        enddo   
!	z12 = (nz-1)*nz*(nz+1)/12.
        do iy = 1,ny
          do ix = 1,nx
            dbdz(ix,iy) = dbdz(ix,iy)/z12
            bavg(ix,iy) = bavg(ix,iy)/ng
          enddo   
        enddo   
        zavg = zavg/ng
        davgdz = davgdz/z12

!	calculate correlation coefficient for each pixel

        sumzsq = 0.
        do iz = 1,nz
          if((wt(iz)==0.).and.(.not.scan)) cycle
          dz = iz-z2
          zi = avgb(iz)-(zavg+dz*davgdz)
          sumzsq = sumzsq+zi**2
          do iy = 1,ny
            do ix = 1,nx
              dbrr = brr(ix,iy,iz)-(bavg(ix,iy)+dz*dbdz(ix,iy))
              corr(ix,iy) = corr(ix,iy)+zi*dbrr
            enddo   
          enddo   
        enddo   
        sumcf = 0.
        sumff = 0.
        do iy = 1,ny
          do ix = 1,nx
            fl = flat(ix,iy)
            if(fl>(0.)) then
              corr(ix,iy) = corr(ix,iy)/sumzsq
              sumcf = sumcf+corr(ix,iy)/fl
              sumff = sumff+1./fl**2
            else
              corr(ix,iy) = 0.
            endif
          enddo   
        enddo   

!	make corr orthogonal to 1/flat

        rat = sumcf/sumff
!	do iy = 1,ny
!         do ix = 1,nx
!	   if(flat(ix,iy)>(0.)) &
!            corr(ix,iy) = corr(ix,iy)-rat/flat(ix,iy)
!        enddo   
!      enddo   

        if(verbose) then
          print '(" Plotting sky noise array")'
          call grey(mx,nx,ny,corr,0.,0.,iwin1)
          call waitasec(pause,.true.,corr,mx,nx,ny)
        endif

!	find amount of corr in each frame and subtract it

!	do iz = 1,nz
!	  if((wt(iz)==0.).and.(.not.scan)) goto 300
!	  dz = iz-z2
!	  suma = 0.
!	  sumb = 0.
!	  sumc = 0.
!	  do iy = 1,ny
!	    do ix = 1,nx
!	      b = bavg(ix,iy)+dz*dbdz(ix,iy)
!	      darr = arr(ix,iy,iz)-b
!	      dbrr = brr(ix,iy,iz)-b
!	      if(verbose) plot(ix,iy) = darr-dbrr
!	      if(flat(ix,iy)==0.) then
!	        c = corr(ix,iy)
!	      else
!	        c = corr(ix,iy)-rat/flat(ix,iy)
!	      endif
!	      suma = suma+c*darr
!	      sumb = sumb+c*dbrr
!	      sumc = sumc+c**2
!            enddo   
!          enddo   
!	  if(verbose) then
!	    call grey(mx,nx,ny,plot,0.,0.,iwin1)
!	    call waitasec(pause,.true.,plot,mx,nx,ny)
!	  endif
!	  rata = suma/sumc
!	  ratb = sumb/sumc
!	  do iy = 1,ny
!	    do ix = 1,nx
!	      arr(ix,iy,iz) = arr(ix,iy,iz)-rata*corr(ix,iy)
!	      brr(ix,iy,iz) = brr(ix,iy,iz)-ratb*corr(ix,iy)
!	      if(verbose) plot(ix,iy) = arr(ix,iy,iz)-brr(ix,iy,iz)
!            enddo   
!          enddo   
!	  if(verbose) then
!	    call grey(mx,nx,ny,plot,0.,0.,iwin1)
!	    call waitasec(pause,.true.,plot,mx,nx,ny)
!	  endif
!        enddo   

!	instead find amount of corr in diff and subtract it from arr

        do iz = 1,nz
          if((wt(iz)==0.).and.(.not.scan)) cycle
          dz = iz-z2
          suma = 0.
          sumc = 0.
          do iy = 1,ny
            do ix = 1,nx
              if((iy>1).and.(iy<ny)) then
                c = (corr(ix,iy-1)+2.*corr(ix,iy)+corr(ix,iy+1))/4.
                d = (arr(ix,iy-1,iz)+2.*arr(ix,iy,iz)+arr(ix,iy+1,iz) &
                 -(brr(ix,iy-1,iz)+2.*brr(ix,iy,iz)+brr(ix,iy+1,iz)))/4.
              else
                c = corr(ix,iy)
                d = arr(ix,iy,iz)-brr(ix,iy,iz)
              endif
              if(verbose) plot(ix,iy) = d
              suma = suma+c*d
              sumc = sumc+c**2
            enddo   
          enddo   
          if(verbose) then
            call grey(mx,nx,ny,plot,0.,0.,iwin1)
            call waitasec(pause,.true.,plot,mx,nx,ny)
          endif
          rat = suma/sumc
          do iy = 1,ny
            do ix = 1,nx
              arr(ix,iy,iz) = arr(ix,iy,iz)-rat*corr(ix,iy)
              if(verbose) plot(ix,iy) = arr(ix,iy,iz)-brr(ix,iy,iz)
            enddo   
          enddo   
          if(verbose) then
            call grey(mx,nx,ny,plot,0.,0.,iwin1)
            call waitasec(pause,.true.,plot,mx,nx,ny)
          endif
        enddo   

        ierr = 0
        return
      end


      subroutine cirrus(arr,brr,flat,nz,ierr)
!	remove sky noise by subtracting (a+by)*brr+(c+dy)/flat from arr
!	with a,b,c,d chosen to minimize the sum of (arr-brr)**2
!	Note: sky noise depends on y because clouds can vary during an
!	array readout.

        use ius
        use dims
        use modes

        real arr(mx,my,1),brr(mx,my,1),flat(mx,my)
        real a(4),x(4),alpha(4,4),beta(4)
!	real plot(mx,my)
        real krot
        logical good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical scan

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /iwin/ iwin1,iwin2

        print '(" Subtracting sky fluctuations in cirrus")'
        print '("*** Warning: cirrus removes average continuum", &
          " along slit")'
        print '(" Set skyint to use subcorr")'
        write(iupl,'(" Subtracting sky fluctuations in cirrus")')
        scan = (modeobs==MSCAN).or.(modeobs==MMAP)

        ni = 4
        xny = ny
        xny2 = (ny+1.)/2.

        do iz = 1,nz

          do j = 1,ni
            do i = 1,ni
              alpha(i,j) = 0.
            enddo   
            beta(j) = 0.
          enddo   
          do iy = 1,ny
            dy = (iy-xny2)/xny
            do ix = 1,nx
              if(flat(ix,iy)<=0.) cycle
              x(1) = brr(ix,iy,iz)
              x(2) = x(1)*dy
              x(3) = 1./flat(ix,iy)
              x(4) = x(3)*dy
              diff = arr(ix,iy,iz)-brr(ix,iy,iz)
              do j = 1,ni
                do i = 1,ni
                  alpha(i,j) = alpha(i,j)+x(i)*x(j)
                enddo   
                beta(j) = beta(j)+diff*x(j)
              enddo   
            enddo   
          enddo   

          call matinv(alpha,ni,det)
          if(det==0.) then
            print '(" Can''t find sky parameters.")'
            write(iupl,'(" Can''t find sky parameters.")')
            ierr = 1
            return
          endif

          do i = 1,ni
            ai = 0.
            do j = 1,ni
              ai = ai+beta(j)*alpha(i,j)
            enddo   
            a(i) = ai
          enddo   

          do iy = 1,ny
            dy = (iy-xny2)/xny
            do ix = 1,nx
              if(flat(ix,iy)<=0.) cycle
              arr(ix,iy,iz) = arr(ix,iy,iz) &
                -(a(1)+a(2)*dy)*brr(ix,iy,iz) &
                -(a(3)+a(4)*dy)/flat(ix,iy)
!	     if(verbose) plot(ix,iy) = arr(ix,iy,iz)-brr(ix,iy,iz)
            enddo   
          enddo   
          if(verbose) then
            print '(" Sky noise parameters: ",4f10.5)', (a(i),i=1,4)
!	   call grey(mx,nx,ny,plot,0.,0.,iwin1)
!	   call waitasec(pause,.true.,plot,mx,nx,ny)
          endif

        enddo   

        ierr = 0
        return
      end


      subroutine setillum(flat,ierr)
!	old setillum based on an untorted flat
!	set illum = -1 if pixel is shifted in from off the array
!	illum = 0 if pixel is between orders
!	or if pixel is shifted from where flat = 0

        use dims
        use modes
        use consts

        real flat(mx,my)
        real nodpa,lores,kmirror,krot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical good

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /iwin/ iwin1,iwin2

        if(modetort==MNONE) then
          do iy = 1,ny
            do ix = 1,nx
              if(flat(ix,iy)==0.) then
                illum(ix,iy) = 0
              else
                illum(ix,iy) = 1
              endif
            enddo   
          enddo   
          return
        endif

!	find which pixels of untorted flat(iu,iv) map onto each pixel of
!	torted illum(ix,iy)

        if(crossdisp) then
          hrslrot = -slitrot+krot
!	include echelon smile and slit rotation here
          hrskew = 2.*hrg*abs(hrr)+tan(hrslrot)
!	include xd smile and spectrum rotation by k mirror here
          xdskew = 2.*xdg*xdr+tan(krot)
          xdsmile = -xdr*pixelwd/xdfl
          xddisp = (xdr*xdfl)/(abs(hrr)*hrfl)
          xdnlin = -(xdr+1./(2.*xdr))*(pixelwd/xdfl)/2.
          hrnlin = -(abs(hrr)+1./(2.*abs(hrr)))*(pixelwd/hrfl)/2.
          xorder0 = xorder1-spacing/2.
          ispac = int(spacing)
          if(ispac/=nt) &
            print '("***ispac, nt in setillum:",2i4)', ispac,nt
          spaci = ispac
        else
          xdskew = 2.*xdg*xdr+tan(slitrot+2.*krot)
          xdsmile = -xdr*pixelwd/xdfl
          xdnlin = -(xdr+1./(2.*xdr))*(pixelwd/xdfl)/2.
        endif
        cosrot = cos(detrot)
        sinrot = sin(detrot)
        xmid = (nx+1)/2.
        ymid = (ny+1)/2.

        do iy = 1,ny
          do ix = 1,nx
            x = ix-xmid
            y = iy-ymid
            if(crossdisp) then
!	slit skewing within orders by echelon smile
!	     order = (ix-xorder0)/spacing
              order = ix/spaci+0.5
              iorder = nint(order)
              dorder = order-iorder
              x = xorder0+order*spacing-xmid
              dx = spacing*dorder
              y = y+hrskew*dx
!	non-linearity of echelon spectrum
              y = y+hrnlin*(y**2-ymid**2)
!	skewing by cross dispersion, K mirror, and cross-dispersion smile
!	note: xd dispersion depends on linear y (wavelength)
              x = x+xddisp*(iy-ymid)+xdskew*y+dxdskew*x*y+xdsmile*y**2
!	non-linearity of xd spectrum
              x = x+xdnlin*(x**2-xmid**2)
            elseif(modeinst/=MCAMERA) then
!	skewing by cross-dispersion smile
              x = x+xdskew*y+xdsmile*y**2
              x = x+xdnlin*(x**2-xmid**2)
            endif
!	array rotation
            u = x*cosrot-y*sinrot
            v = y*cosrot+x*sinrot
!	check status of 4 neighboring pixels
            iu1 = int(u+xmid)
            iu2 = iu1+1
            iv1 = int(v+ymid)
            iv2 = iv1+1
            if((iu1>1).and.(iu2<nx) &
              .and.(iv1>1).and.(iv2<ny)) then
              if((modeinst==MMED).or.(modeinst==MLOW)) then
                illum(ix,iy) = illx(iy)
              else
                illum(ix,iy) = illx(ix)
              endif
              if((flat(iu1,iv1)==0.).or.(flat(iu2,iv1)==0.) &
                .or.(flat(iu1,iv2)==0.).or.(flat(iu2,iv2)==0.)) &
                illum(ix,iy) = 0
            else
              illum(ix,iy) = -1
            endif
          enddo   
        enddo

        return
      end


      subroutine midillum(blk,ierr)
!	new setillum based on a torted blk or blk-sky
!	set illum = -1 if pixel is shifted in from off the array
!	illum = 0 if pixel is between orders
!	or if pixel is shifted from where flat = 0
!	also set illx = 1 if >ny/2 pixels are illuminated in a column

        use dims
        use modes
        use consts

        real blk(mx,my),pow(mx)
        real nodpa,lores,kmirror,krot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical good

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /iwin/ iwin1,iwin2

        if(crossdisp) then
          powmin = 1.0e+12
          powmax = 0.
          ix = 0
          do it = 1,nt
            nill = 0
            powi = 0.
            do iorder = 1,norder
              ix = (iorder-1)*nt+it
              do iy = 1,ny
                if(blk(ix,iy)>0.) then
                  nill = nill+1
                  powi = powi+blk(ix,iy)
                endif
              enddo
            enddo
            if(nill>norder*ny/2) then
              powi = powi/nill
              do iorder = 1,norder
                ix = (iorder-1)*nt+it
                pow(ix) = powi
              enddo
              if(powi<powmin) powmin = powi
              if(powi>powmax) powmax = powi
            else
              do iorder = 1,norder
                ix = (iorder-1)*nt+it
                pow(ix) = 0.
              enddo
            endif
          enddo
          do jx = ix+1,nx
            pow(jx) = 0.
          enddo
          thresh = powmin+thrfac*(powmax-powmin)
          do ix = 1,nx
            if(pow(ix)==0.) then
              illx(ix) = -1
            elseif(pow(ix)<thresh) then
              illx(ix) = 0
            else
              illx(ix) = 1
            endif
          enddo
        else
          do iy = 1,ny
            nill = 0
            powi = 0.
            do ix = 1,nx
              if(blk(ix,iy)>0.) then
                nill = nill+1
                powi = powi+blk(ix,iy)
              endif
            enddo
            if(nill>nx/2) then
              powi = powi/nill
              pow(iy) = powi
              if(powi<powmin) powmin = powi
              if(powi>powmax) powmax = powi
            else
              pow(iy) = 0.
            endif
          enddo
          thresh = powmin+thrfac*(powmax-powmin)
          do iy = 1,ny
            if(pow(iy)==0.) then
              illx(iy) = -1
            elseif(pow(iy)<thresh) then
              illx(iy) = 0
            else
              illx(iy) = 1
            endif
          enddo
        endif

        do iy = 1,ny
          kx = 0
          lx = 0
          do ix = 1,nx
            if(blk(ix,iy)/=0.) then
              illum(ix,iy) = 1
              if(kx==0) kx = ix
              lx = ix
            else
              if(kx==0) then
!	haven''t yet seen an illuminated pixel in this row
                illum(ix,iy) = -1
              else
                illum(ix,iy) = 0
              endif
            endif
          enddo
!	pixels beyond last illuminated pixel
          do ix = lx+1,nx
            illum(ix,iy) = -1
          enddo
        enddo
        do ix = 1,nx
          ky = 0
          ly = 0
          illsum = 0
          do iy = 1,ny
            if(illum(ix,iy)>0) then
              illsum = illsum+1
              if(ky==0) ky = iy
              ly = iy
            endif
          enddo
          do iy = 1,ky-1
            illum(ix,iy) = -1
          enddo
          do iy = ly+1,ny
            illum(ix,iy) = -1
          enddo
        enddo

        return
      end


      subroutine newillum(blk,ierr)
!	new setillum based on a torted blk or blk-sky
!	set illum = -1 if pixel is shifted in from off the array
!	illum = 0 if pixel is between orders in illx
!	or if pixel is shifted from where flat = 0

        use dims
        use modes
        use consts

        real blk(mx,my),pow(mx)
        real nodpa,lores,kmirror,krot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical good

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /iwin/ iwin1,iwin2

        do iy = 1,ny
          kx = 0
          lx = 0
          do ix = 1,nx
            if(blk(ix,iy)/=0.) then
              illum(ix,iy) = 1
              if(kx==0) kx = ix
              lx = ix
            else
              if(kx==0) then
!	haven''t yet seen an illuminated pixel in this row
                illum(ix,iy) = -1
              else
                illum(ix,iy) = 0
              endif
            endif
          enddo
!	pixels beyond last illuminated pixel
          do ix = lx+1,nx
            illum(ix,iy) = -1
          enddo
        enddo
        do ix = 1,nx
          ky = 0
          ly = ny
          do iy = 1,ny
            if(illum(ix,iy)>0) then
              if(ky==0) ky = iy
              ly = iy
            endif
          enddo
          do iy = 1,ky-1
            illum(ix,iy) = -1
          enddo
          do iy = ly+1,ny
            illum(ix,iy) = -1
          enddo
        enddo

        if(crossdisp) then
          do ix = 1,nx
            if(illx(ix)==0) then
              do iy = 1,ny
                illum(ix,iy) = 0
              enddo
            endif
          enddo
        else
          do iy = 1,ny
            if(illx(iy)==0) then
              do ix = 1,nx
                illum(ix,iy) = 0
              enddo
            endif
          enddo
        endif

        return
      end


      subroutine maketempl(arr,flat,std,templ,nz, &
          submean,flatnorm,sumtempl,ierr)
!	make a template array for shift, wtadd, and extract
!	sum arr in z and average spectrally, weighted by 1/flat**2 or 1/std**2
!	if flatnorm = true, divide spectral sum by sum of weights
!	else divide spectral sum by average weight sum
!	so different spatial points retain different weighting

        use ius
        use dims
        use modes

        real arr(mx,my,1),flat(mx,my),std(mx,my),templ(mx,my)
        logical submean,flatnorm,sumtempl
        real templt(mx),templx(mx),tmplmed(mx)
        real nodpa,lores,kmirror,krot
        logical good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical tspec,tarr
        character(8) yn
        character(60) comment

        common /goods/good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime
        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /iwin/ iwin1,iwin2

!	average arr in z

        tspec = sumtempl.and.crossdisp
        tarr = sumtempl.and.(.not.crossdisp)
        if(nz==1) wt(1) = 1.
        wtsum = 0.
        do iz = 1,nz
          if(wt(iz)/=0.) wtsum = wtsum+1.
        enddo   
        if(wtsum==0.) then
          print '("*** All weights = 0.  Can''t make template")'
          write(iupl,'("*** All weights = 0.  Can''t make template")')
          ierr = 9
          return
        endif
        nbad = 0
        do iy = 1,ny
          do ix = 1,nx
            sum = 0.
            if((.not.tarr).and.(illum(ix,iy)==1)) then
              do iz = 1,nz
                sum = sum+wt(iz)*arr(ix,iy,iz)
              enddo   
            endif
            if(tarr.and.(iy>4)) then
              templ(ix,iy) = arrsum(ix,iy)
            elseif(sum==sum) then
              templ(ix,iy) = sum/wtsum
            else
              if(nbad<4) print'(" templ(",2i4,") =",e12.3)',ix,iy,sum
              nbad = nbad+1
              templ(ix,iy) = 0.
            endif
          enddo   
        enddo   
        if(nbad>4) print '(i8," bad points in template")', nbad

!	collapse and replicate template spectrally (over y and orders)
        if((intshift(3)*intshift(4))>0) then
          iorder1 = max(1,intshift(3))
          iorder2 = min(norder,intshift(4))
        else
          iorder1 = 1
          iorder2 = norder
        endif
        if((intshift(5)*intshift(6))>0) then
          if(crossdisp) then
            iy1 = max(1,intshift(5))
            iy2 = min(ny,intshift(6))
          else
            ix1 = max(1,intshift(5))
            ix2 = min(nx,intshift(6))
          endif
        else
          if(crossdisp) then
            iy1 = 1
            iy2 = ny
          else
            ix1 = 1
            ix2 = nx
          endif
        endif

        itext1 = intext(1)
        itext2 = intext(2)
        itext3 = intext(3)
        itext4 = intext(4)
        if(modeext<MINTWT) then
          itext1 = 0
          itext2 = 0
        elseif(itext3*itext4>0) then
          if(itext3<itext1) itext1 = itext3
          if(itext4>itext2) itext2 = itext4
        endif

        if(crossdisp) then

          ispac = int(spacing)
          if(ispac/=nt) print '("***ispac, nt in maketempl:",2i4)', &
            ispac,nt

!	collapse templ or specsum into templt
          sumsum = 0.
          illsumsum = 0
          if((intext(1)>nt).or.(intext(1)>intext(2))) then
            print '("*** All pixels outside of extint.", &
              "  Hit return to continue",$)'
            read '(a8)', yn
            ierr = 9
            return
          endif
          do it = 1,nt
            if((itext1*itext2>0) &
              .and.((it<itext1).or.(it>itext2))) then
              templt(it) = 0.
              cycle
            endif
            its = it+4
            if(it==nt) its = nt+3
            sum = 0.
            swtsum = 0.
            illsum = 0
            do iorder = iorder1,iorder2
!	     ix0 = int(xorder1+(iorder-1)*spacing)
              ix0 = (iorder-1)*nt
              is0 = ny*(iorder-1)
              ix = ix0+it
              if(ix>nx) exit
              do iy = iy1,iy2
                if(illum(ix,iy)==1) then
!	          if(stdwt) then
!	            swt = 1./std(ix,iy)**2
!	          else
                    swt = 1./flat(ix,iy)**2
!	          endif
                  if(tspec) then
                    is = is0+iy
                    sum = sum+swt*specsum(is,its)
                  else
                    sum = sum+swt*templ(ix,iy)
                  endif
                  swtsum = swtsum+swt
                  illsum = illsum+1
                endif
              enddo   
            enddo   
            if(illsum>((iorder2-iorder1+1)*(iy2-iy1+1)/2)) then
              if(flatnorm) then
                sum = sum/swtsum
                sumsum = sumsum+sum
              else
                sumsum = sumsum+swtsum
              endif
              illsumsum = illsumsum+1
            else
              sum = 0.
            endif
            templt(it) = sum
          enddo   

          if((sumsum==0.).or.(illsumsum==0)) then
            print '("*** All values = 0 in template")'
            print '(" Hit RETURN to continue",$)'
            read '(a8)', yn
            ierr = 9
            return
          endif
          if(flatnorm) then
            avg = sumsum/illsumsum
          else
            swtsum = sumsum/illsumsum
            sumsum = 0.
            illsum = 0
            do it = 1,nt
              if(templt(it)/=0.) then
                templt(it) = templt(it)/swtsum
                sumsum = sumsum+templt(it)
                illsum = illsum+1
              endif
            enddo   
            avg = sumsum/illsum
          endif

!	subtract mean for use with addwt and extract if dosubsky
          if(submean) then
!	try subtracting median instead
            if(nodon) then
              call median(tmplmed,templt,1,1,nt,2,ierr)
            else
              call median(tmplmed,templt,1,1,nt,1,ierr)
            endif
            do it = 1,nt
              if(templt(it)/=0.) templt(it) = templt(it)-tmplmed(1)
            enddo   
          endif

!	print templt for extraction weighting
          if(nz==1) then
            print '(" Extraction template:")'
            print '(8f8.3)', (templt(it), it=1,nt)
            write(iupl,'(" Extraction template:")')
            write(iupl,'(8f8.3)') (templt(it), it=1,nt)
            if(irtf.and.nodon.and.(slitpa>720.)) then
              itmin = 0
              itmax = 0
              templmin = 0.
              templmax = 0.
              do it = 1,nt
                templi = templt(it)
                if(templi<templmin) then
                  templmin = templi
                  itmin = it
                elseif(templi>templmax) then
                  templmax = templi
                  itmax = it
                endif
              enddo
              if(itmin<itmax) then
                slitpa = nodpa
              elseif(itmin>itmax) then
                slitpa = mod(nodpa+180.,360.)
              else
                print '(" All template values = 0.")'
              endif
              print '(" Setting slitpa =",f8.2)', slitpa
              write(iupl,'(" Setting slitpa =",f8.2)') slitpa
              comment = 'reset based on AB orientation and nod'
              if(dofits) call fithreal('SLITPA  ',slitpa,comment,iufith)
              write(iuredh,'("slitpa  = ",f8.2)') slitpa
              if(dosum) write(iusumh,'("slitpa  = ",f8.2)') slitpa
            endif
          endif

!	find edges of illumination of orders (ends of slit)
          templ1 = 0.
          it1 = 0
          do it = 1,nt
            if(templt(it)/=0.) then
              if(it1==0) then
                templ1 = templt(it)
                it1 = it
              endif
              templ2 = templt(it)
              it2 = it
            endif
          enddo   

!	replicate templt over orders
!	include end orders if only unilluminated pixels fall off array
            do ix = 1,nx
!	     iorder = (ix-xorder1)/spacing+1.
              iorder = (ix-1)/nt+1
              it = ix-(iorder-1)*nt
              if((iorder<1).or.(iorder>norder).or.(illx(ix)<1)) then
                templx(ix) = 0.
              else
                templx(ix) = templt(it)
              endif
          enddo   

!	replicate templx over y
          do iy = 1,ny
            do ix = 1,nx
              if(illum(ix,iy)>=0) then
                templ(ix,iy) = templx(ix)
              else
                templ(ix,iy) = 0.
              endif
            enddo   
          enddo   

        elseif((modeinst==MMED).or.(modeinst==MLOW)) then

          if((itext1*itext2>0).and.(itext2<=ny)) then
            iy1 = itext1
            iy2 = min(ny,itext2)
          else
            iy1 = 1
            iy2 = ny
          endif
          sumsum = 0.
          illrow = 0
          do iy = 1,ny
            templt(iy) = 0.
          enddo   
          do iy = iy1,iy2
            sum = 0.
            swtsum = 0.
            illsum = 0
            do ix = ix1,ix2
              if(illum(ix,iy)==1) then
!	        if(stdwt) then
!	          swt = 1./std(ix,iy)**2
!	        else
                  swt = 1./flat(ix,iy)**2
!	        endif
                sum = sum+swt*templ(ix,iy)
                swtsum = swtsum+swt
                illsum = illsum+1
              endif
            enddo   
            if(illsum>((ix2-ix1+1)/2)) then
              if(flatnorm) then
                sum = sum/swtsum
                sumsum = sumsum+sum
              else
                sumsum = sumsum+swtsum
              endif
              illrow = illrow+1
            else
              sum = 0.
            endif
            templt(iy) = sum
          enddo   
          if(irtf.and.nodon.and.(slitpa>720.)) then
            itmin = 0
            itmax = 0
            templmin = 0.
            templmax = 0.
            do it = 1,ny
              templi = templt(it)
              if(templi<templmin) then
                templmin = templi
                itmin = it
              elseif(templi>templmax) then
                templmax = templi
                itmax = it
              endif
            enddo
            if(itmin<itmax) then
              slitpa = nodpa
            elseif(itmin>itmax) then
              slitpa = mod(nodpa+180.,360.)
            else
              print '(" All template values = 0.")'
            endif
            print '(" Setting slitpa =",f6.2)', slitpa
            write(iupl,'(" Setting slitpa =",f6.2)') slitpa
            if(dofits) call fithreal('SLITPA  ',slitpa,comment,iufith)
            write(iuredh,'("slitpa  = ",f12.4)') slitpa
            if(dosum) write(iusumh,'("slitpa  = ",f12.4)') slitpa
          endif

          if((sumsum==0.).or.(illrow==0)) then
            print '("*** All values = 0 in template")'
            print '(" Hit RETURN to continue",$)'
            read '(a8)', yn
            ierr = 9
            return
          endif
          if(flatnorm) then
            avg = sumsum/illrow
          else
            swtsum = sumsum/illrow
            sumsum = 0.
            illsum = 0
            do iy = iy1,iy2
              if(templt(iy)/=0.) then
                templt(iy) = templt(iy)/swtsum
                sumsum = sumsum+templt(iy)
                illsum = illsum+1
              endif
            enddo   
            avg = sumsum/illsum
          endif

          if(submean) then
            do iy = iy1,iy2
              if(templt(iy)/=0.) templt(iy) = templt(iy)-avg
            enddo   
          endif

          do iy = 1,ny
            do ix = 1,nx
              if(illum(ix,iy)==1) then
                templ(ix,iy) = templt(iy)
              else
                templ(ix,iy) = 0.
              endif
            enddo   
          enddo   

        endif

        if(verbose) then
          print '(" Plotting template")'
          call grey(mx,nx,ny,templ,0.,0.,iwin1)
          call waitasec(pause,.true.,templ,mx,nx,ny)
        endif

        ierr = 0
        return
      end


      subroutine shift(arr,flat,std,nz,ierr)
!	correlate array of frames with sum and shift frames spatially
!	if(dosum) correlate with sumspec
!	if(intshift(1)<0) apply intshift(2) instead of finding shift

        use ius
        use dims
        use modes

        real arr(mx,my,mp),flat(mx,my),std(mx,my)
        real templ(mx,my),acoll(mx),tcoll(mx),scoll(mx),corr(17)
        real sumarr(mx,my),plot(mx,my)
        real nodpa,lores,kmirror,krot
        integer isha(mp)
        logical good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical submean,flatnorm,sumtempl,fixshft
        logical hdopen
        character(8) yn
        character(60) comment

        common /nn/ nx,ny,nc,norder,ns,nt
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /iwin/ iwin1,iwin2

        if((.not.doshift).or.(nz<2)) then
          print '(" Can''t shift one frame")'
          return
        endif

        sumwt = 0.
        do iz = 1,nz
          sumwt = sumwt+abs(wt(iz))
        enddo   
        if(sumwt==0.) then
          print '("!!! All weights are zero in shift")'
          write(iupl,'("!!! All weights are zero in shift")')
          return
        endif
        xnz = sumwt

        if(intshift(1)<0) then
          fixshft = .true.
          ishft = intshift(2)
          print '(" Shifting all pairs by ",i3)', ishft
          do iz = 1,nz
            isha(iz) = ishft
          enddo   
          goto 50
        else
          fixshft = .false.
        endif

        sumsum = 0.
        do iy = 1,ny
          do ix = 1,nx
            sum = 0.
            if(illum(ix,iy)==1) then
              do iz = 1,nz
                sum = sum+wt(iz)*arr(ix,iy,iz)
              enddo   
            endif
            sumarr(ix,iy) = sum/sumwt
            sumsum = sumsum+sum/sumwt
          enddo   
        enddo   

!	make template from sum of frames
!	2/02 don't require dosubsky = true for submean = true

        if((modeinst==MHIMED).or.(modeinst==MHILOW)) then
          slitlen = 0.8*spacing
        elseif((modeinst==MMED).or.(modeinst==MLOW)) then
          slitlen = 100.
        else
          slitlen = 0.
        endif
        submean = (((modeobs==MNOD).or.(modeobs==MMAP)) &
          .and.(dist<slitlen))
        flatnorm = .false.
        sumtempl = (dosum.and.(sumtime>0))
        call maketempl(arr,flat,std,templ,nz, &
          submean,flatnorm,sumtempl,ierr)

!	cross-dispersed or camera mode; shift in x

   50 if((modeinst==MHIMED).or.(modeinst==MHILOW) &
        .or.(modeinst==MCAMERA)) then
        if(fixshft) goto 250
        nsh = min0(4,nt/3)
        ncorr = 2*nsh+1


!	collapse template and stddev in y, weighting with 1/flat**2

        if((intshift(1)>0).and.(intshift(2)>=intshift(1))) then
          it1 = max(1,intshift(1))
          it2 = min(nt,intshift(2))
        else
          it1 = 1
          it2 = nt
        endif
        if((intshift(3)>0).and.(intshift(4)>=intshift(3))) then
          iorder1 = max(1,intshift(3))
          iorder2 = min(norder,intshift(4))
        else
          iorder1 = 1
          iorder2 = norder
        endif
        if((intshift(5)>0).and.(intshift(6)>=intshift(5))) then
          iy1 = max(1,intshift(5))
          iy2 = min(ny,intshift(6))
        else
          iy1 = 1
          iy2 = ny
        endif
        do ix = 1,nx
          sum = 0.
          sumstd = 0.
          sumwt = 0.
          it = mod(ix,nt)
          iorder = (ix-1)/nt+1
          if((iorder<iorder1).or.(iorder>iorder2) &
            .or.(it<it1).or.(it>it2)) goto 210
          do iy = iy1,iy2
            if(illum(ix,iy)>0) then
              flwt = 1./flat(ix,iy)**2
              sum = sum+templ(ix,iy)*flwt
              sumstd = sumstd+std(ix,iy)*flwt
              sumwt = sumwt+flwt
            endif
          enddo   
  210     if(sumwt>0.) then
            tcoll(ix) = sum/sumwt
            scoll(ix) = sumstd/sumwt
          else
            tcoll(ix) = 0.
            scoll(ix) = 1.
          endif
          plot(ix,1) = tcoll(ix)
        enddo   

        msh = 0
        nmsh = 0
        do iz = 1,nz

!	collapse data array in y, weighting with 1/flat**2

          do ix = 1,nx
            sum = 0.
            sumwt = 0.
            it = mod(ix,nt)
            iorder = (ix-1)/nt+1
            if((iorder<iorder1).or.(iorder>iorder2) &
              .or.(it<it1).or.(it>it2)) goto 220
            do iy = iy1,iy2
              if(illum(ix,iy)>0) then
                flwt = 1./flat(ix,iy)**2
                sum = sum+arr(ix,iy,iz)*flwt
                sumwt = sumwt+flwt
              endif
            enddo   
  220       if(sumwt>0.) then
              sum = sum/sumwt
            else
              sum = 0.
            endif
            plot(ix,iz+1) = sum
            acoll(ix) = sum*wt(iz)
          enddo   
          if(wt(iz)==0.) cycle

!	find shift which maximizes contribution to S**2 or(S/N)**2

          do ish = -nsh,nsh
            icorr = ish+nsh+1
            sum = 0.
            do ix = 1+nsh,nx-nsh
              ixsh = ix+ish
              sum = sum+(xnz*tcoll(ix)+wt(iz)*acoll(ixsh))**2 &
                /(xnz*scoll(ix)**2+(wt(iz)*scoll(ixsh))**2)
            enddo   
            corr(icorr) = sum
          enddo   

          imax = 1
          cmax = corr(1)
          do icorr = 2,ncorr
            if(corr(icorr)>cmax) then
              cmax = corr(icorr)
              imax = icorr
            endif
          enddo   
          icorr = imax
          ish = icorr-nsh-1
          isha(iz) = ish
          if((icorr==1).or.(icorr==ncorr)) then
            wt(iz) = 0.
          else
            msh = msh+ish
            nmsh = nmsh+1
          endif
        enddo   

        if(nmsh==0) then
          print '("***All x shifts are out of range")'
          if(dosum) print '(" Try starting a new sumspec")'
          print '(" Shift parameters:")'
          print '(8i8)', (-isha(i),i=1,nz)
          print '(8es10.3)', (corr(icorr),icorr=1,2*nsh+1)
          print '(" Hit RETURN to continue",$)'
          read '(a8)', yn
          write(iupl,'("***All x shifts are out of range")')
          write(iupl,'(" Shift parameters:")')
          write(iupl,'(8i8)') (-isha(i),i=1,nz)
          ierr = 1
          goto 900
        endif
        msh = msh/nmsh
        if((.not.dosum).and.(iabs(msh)>0)) then
          do iz = 1,nz
            if(wt(iz)/=0.) isha(iz) = isha(iz)-msh
          enddo   
        endif
        print '(" Shift parameters:")'
        print '(8i8)', (-isha(i),i=1,nz)
        write(iupl,'(" Shift parameters:")')
        write(iupl,'(8i8)') (-isha(i),i=1,nz)

!	shift array spatially

  250   do iz = 1,nz
          ish = isha(iz)
          if((wt(iz)==0.).or.(ish==0)) then
            continue
          elseif(ish>0) then
            ix1 = 1
            ix2 = nx-ish
            do ix = ix1,ix2
              ixsh = ix+ish
              do iy = 1,ny
                if(illum(ix,iy)==1) then
                  if(illum(ixsh,iy)==1) then
                    arr(ix,iy,iz) = arr(ixsh,iy,iz)
                  else
                    arr(ix,iy,iz) = sumarr(ix,iy)
                  endif
                else
                  arr(ix,iy,iz) = 0.
                endif
              enddo   
            enddo   
          elseif(ish<0) then
            ix1 = nx
            ix2 = 1-ish
            do ix = ix1,ix2,-1
              ixsh = ix+ish
              do iy = 1,ny
                if(illum(ix,iy)==1) then
                  if(illum(ixsh,iy)==1) then
                    arr(ix,iy,iz) = arr(ixsh,iy,iz)
                  else
                    arr(ix,iy,iz) = sumarr(ix,iy)
                  endif
                else
                  arr(ix,iy,iz) = 0.
                endif
              enddo   
            enddo   
          endif
        enddo   

      endif

!	long-slit or camera mode; shift in y

      if((modeinst==MMED).or.(modeinst==MLOW) &
        .or.(modeinst==MCAMERA)) then
        if(fixshft) goto 150
        nsh = 4
        ncorr = 2*nsh+1

        iy1 = intshift(1)
        iy2 = intshift(2)
        if((iy1<1).or.(iy1>ny)) then
          iy1 = 1
          iy2 = ny
        elseif((iy2<iy1).or.(iy2>ny)) then
          iy2 = ny
        endif
        if((intshift(5)*intshift(6))>0) then
          ix1 = max(1,intshift(5))
          ix2 = min(nx,intshift(6))
        else
          ix1 = 1
          ix2 = nx
        endif

!	collapse template and stddev spectrally, weighting with 1/flat**2

        do iy = 1,ny
          sum = 0.
          sumstd = 0.
          sumflat = 0.
          do ix = ix1,ix2
            if(illum(ix,iy)==1) then
              flwt = 1./flat(ix,iy)**2
              sum = sum+templ(ix,iy)*flwt
              sumstd = sumstd+std(ix,iy)*flwt
              sumflat = sumflat+flwt
            endif
          enddo   
          if(sumflat>0.) then
            tcoll(iy) = sum/sumflat
            scoll(iy) = sumstd/sumflat
          else
            tcoll(iy) = 0.
            scoll(iy) = 1.
          endif
          plot(iy,1) = tcoll(iy)
        enddo   

        msh = 0
        nmsh = 0
        do iz = 1,nz
          if(wt(iz)==0.) then
            do iy = 1,ny
              plot(iy,iz+1) = 0.
            enddo
            cycle
          endif

!	collapse data array spectrally, weighting with 1/flat

          do iy = iy1,iy2
            sum = 0.
            sumflat = 0.
            do ix = ix1,ix2
              if(illum(ix,iy)==1) then
                flwt = 1./flat(ix,iy)
                sum = sum+arr(ix,iy,iz)*flwt
                sumflat = sumflat+flwt
              endif
            enddo   
            if(sumflat>0.) then
              acoll(iy) = wt(iz)*sum/sumflat
            else
              acoll(iy) = 0.
            endif
            plot(iy,iz+1) = acoll(iy)
          enddo   

!	calculate S**2 or (S/N)**2 for different shifts

          do ish = -nsh,nsh
            icorr = ish+nsh+1
            sum = 0.
            do iy = iy1+nsh,iy2-nsh
              iysh = iy+ish
              sum = sum+(tcoll(iy)+acoll(iysh))**2 &
                /(xnz*scoll(iy)**2+(wt(iz)*scoll(iysh))**2)
            enddo   
            corr(icorr) = sum
          enddo   

!	choose shift which maximizes S or S/N

          imax = 1
          cmax = corr(1)
          do icorr = 2,ncorr
            if(corr(icorr)>cmax) then
              cmax = corr(icorr)
              imax = icorr
            endif
          enddo   
          icorr = imax
          ish = icorr-nsh-1
          isha(iz) = ish
          if((icorr==1).or.(icorr==ncorr)) then
            wt(iz) = 0.
          else
            msh = msh+ish
            nmsh = nmsh+1
          endif
        enddo   

        if(nmsh==0) then
          print '("***All y shifts are out of range")'
          print '(" Shift parameters:")'
          print '(8i8)', (-isha(i),i=1,nz)
          print '(" Hit RETURN to continue",$)'
          read '(a8)', yn
          write(iupl,'("***All y shifts are out of range")')
          write(iupl,'(" Shift parameters:")')
          write(iupl,'(8i8)') (-isha(i),i=1,nz)
          ierr = 1
          goto 900
        endif
        print '(" Shift parameters:")'
        print '(8i8)', (-isha(i),i=1,nz)
        write(iupl,'(" Shift parameters:")')
        write(iupl,'(8i8)') (-isha(i),i=1,nz)

!	shift array spatially

  150   do iz = 1,nz
          ish = isha(iz)
          if((wt(iz)==0.).or.(ish==0)) then
            cycle
          elseif(ish>0) then
            do iy = iy1,iy2-ish
              iysh = iy+ish
              do ix = 1,nx
                if(illum(ix,iy)==1) then
                  if(illum(ix,iysh)==1) then
                    arr(ix,iy,iz) = arr(ix,iysh,iz)
                  else
                    arr(ix,iy,iz) = sumarr(ix,iy)
                  endif
                else
                  arr(ix,iy,iz) = 0.
                endif
              enddo   
            enddo   
          elseif(ish<0) then
            do iy = iy2,iy1-ish,-1
              iysh = iy+ish
              do ix = 1,nx
                if(illum(ix,iy)==1) then
                  if(illum(ix,iysh)==1) then
                    arr(ix,iy,iz) = arr(ix,iysh,iz)
                  else
                    arr(ix,iy,iz) = sumarr(ix,iy)
                  endif
                else
                  arr(ix,iy,iz) = 0.
                endif
              enddo   
            enddo   
          endif
        enddo   

      endif

        ierr = 0
        if(fixshft) goto 910
  900   nz1 = nz+1
        scale = -1.
        if(ask) scale = 0.
        print '(" Plotting spatial distributions before shift")'
        call perspec(nx,nz1,plot,scale,iwin2)
        call waitasec(pause,.false.,plot,mx,mx,mx)
  910   inquire(unit=iuredh,opened=hdopen)
        if(hdopen) then
          write(iuredh,'("shift   = T")')
          comment = 'data have been shifted along slit'
          if(dofits) call fithlog('SHIFT   ',.true.,comment,iufith)
          if(dosum.and.(sumtime==0.)) then
            write(iusumh,'("shift   = T")')
            if(dofits) call fithlog('SHIFT   ',.true.,comment,iufish)
          endif
        endif

        return
      end


      subroutine wtadd(arr,sumarr,flat,std,nz,ierr)
!	weight diffs, from noise-weighted correlation with straight sum,
!	and add

        use ius
        use dims
        use modes

        real arr(mx,my,mp),sumarr(mx,my),flat(mx,my),std(mx,my)
        real templ(mx,my),wt(mp),plot(mx,my)
        real nodpa,lores,kmirror,krot
        logical good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical submean,flatnorm,sumtempl
        logical hdopen
        character(60) comment

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt0(mp)
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /iwin/ iwin1,iwin2

        ierr = 0

        if((.not.doaddwt).or.(nz<4)) then
          print '(" Doing unweighted addition of pairs")'
          write(iupl,'(" Doing unweighted addition of pairs")')
          call adddiffs(arr,sumarr,std,nz,ierr)
          addtime = beamtime*nz
          return
        endif

        xnz = 0.
        do iz = 1,nz
          if(wt0(iz)/=0.) xnz = xnz+1.
        enddo   
        if(xnz==0.) then
          print '(" All weights = 0.  Can''t add pairs")'
          write(iupl,'(" All weights = 0.  Can''t add pairs")')
          ierr = 1
          return
        endif

!	make template from straight sum of frames,
!	collapsed and replicated spectrally

        submean = .true.
        flatnorm = .false.
        sumtempl = .false.
        call maketempl(arr,flat,std,templ,nz, &
          submean,flatnorm,sumtempl,ierr)
        if(ierr/=0) then
          print '("*** Error in template.  " &
                "Doing unweighted addition of pairs")'
          write(iupl,'("!!! Error in template.  " &
                "Doing unweighted addition of pairs")')
          call adddiffs(arr,sumarr,std,nz,ierr)
          addtime = beamtime*nz
          ierr = 2
          return
        endif

!	weight in proportion to S/N in spectrum extracted by multiplying
!	by template

        sumwt = 0.
        wtmax = 0.
        do ix = 1,nx
          plot(ix,1) = 0.
        enddo   
        do iz = 1,nz

          if(wt0(iz)==0.) then
            wt(iz) = 0.
            do ix = 1,nx
              plot(ix,iz+1) = 0.
            enddo   
            cycle
          endif

          sum1 = 0.
          sum2 = 0.
          wtsgn = sign(1.,wt0(iz))
          if((modeinst==MHILOW).or.(modeinst==MHIMED)) then
            do ix = 1,nx
              sum3 = 0.
              do iy = 1,ny
                if(illum(ix,iy)==1) then
                  t = templ(ix,iy)
                  a = arr(ix,iy,iz)*wtsgn
!	          if(stdwt) then
!	            s2 = std(ix,iy)**2
!	          else
                    s2 = flat(ix,iy)**2
!	          endif
                  if(s2==0.) then
                    print '("*** Error in wtadd")'
                    print '(" std(",2i4,") =",e10.3)',ix,iy,std(ix,iy)
                    print '(" flat =",e10.3)', flat(ix,iy)
                    cycle
                  endif
!	note: a is shifted, but s2 is not
                  sum1 = sum1+(t*a)/s2
                  sum2 = sum2+(t**2)/s2
                  sum3 = sum3+a/s2
                endif
              enddo   
              plot(ix,1) = plot(ix,1)+sum3/xnz
              plot(ix,iz+1) = sum3
            enddo   
          else
            do iy = 1,ny
              sum3 = 0.
              do ix = 1,nx
                if(illum(ix,iy)==1) then
                  t = templ(ix,iy)
                  a = arr(ix,iy,iz)*wtsgn
!	          if(stdwt) then
!	            s2 = std(ix,iy)**2
!	          else
                    s2 = flat(ix,iy)**2
!	          endif
                  if(s2==0.) then
                    print '("*** Error in wtadd")'
                    print '(" std(",2i4,") =",e10.3)',ix,iy,std(ix,iy)
                    print '(" flat =",e10.3)', flat(ix,iy)
                    cycle
                  endif
                  sum1 = sum1+(t*a)/s2
                  sum2 = sum2+(t**2)/s2
                  sum3 = sum3+a/s2
                endif
              enddo   
              plot(iy,1) = plot(iy,1)+sum3/xnz
              plot(iy,iz+1) = sum3
            enddo   
          endif
          wti = sum1/sum2
          if(wti==0.) then
            print '(" Correlation zero on pair",i3)', iz
            write(iupl,'(" Correlation zero on pair",i3)') iz
            wti = 0.
          elseif(wti<0.) then
            print '(" Correlation negative on pair",i3)', iz
            write(iupl,'(" Correlation negative on pair",i3)') iz
            wti = 0.
          elseif(wti<(0.10)) then
            print '(" Correlation weak on pair",i3,f10.3)', &
              iz,wti
            write(iupl,'(" Correlation weak on pair",i3)') iz
            wti = 0.
          endif
          wt(iz) = wti*wtsgn
          sumwt = sumwt+wti
          if(wti>wtmax) wtmax = wti

        enddo   

        nzero = 0
        do iz = 1,nz
          if(wt(iz)==0.) nzero = nzero+1
          wt(iz) = wt(iz)/sumwt
        enddo   
        wtmax = wtmax/sumwt

        if(nzero>0) then
          print '(i3," pairs given zero weight")', nzero
          write(iupl,'(i3," pairs given zero weight")') nzero
        endif
        if(nzero==nz) then
          print '("*** All weights = 0.  Can''t weight pairs")'
          write(iupl,'("!!! All weights = 0.  Can''t weight pairs")')
          call adddiffs(arr,sumarr,std,nz,ierr)
          addtime = 0.
          ierr = 2
          return
        endif

!	multiply data array by weights and sum over z

        sumwt = 0.
        sumwtsq = 0.
        nzero = 0
        do iz = 1,nz
          if(wt(iz)/=0.) then
            sumwt = sumwt+abs(wt(iz))
            sumwtsq = sumwtsq+wt(iz)**2
          else
            nzero = nzero+1
          endif
        enddo   

        stdrat = sqrt(sumwtsq)/sumwt
        do iy=1,ny
          do ix=1,nx
            if(illum(ix,iy)==1) then
              suma = 0.
              do iz=1,nz
                if(wt(iz)/=0.) then
                  suma = suma+wt(iz)*arr(ix,iy,iz)
                endif
              enddo   
              sumarr(ix,iy) = suma/sumwt
              std(ix,iy) = std(ix,iy)*stdrat
            else
              sumarr(ix,iy) = 0.
              std(ix,iy) = 0.
            endif
          enddo   
        enddo   

        print '(" Weights:")'
        print '(8f8.3)', (wt(i),i=1,nz)
        write(iupl,'(" Weights:")')
        write(iupl,'(8f8.3)') (wt(i),i=1,nz)
        inquire(unit=iuredh,opened=hdopen)
        if(hdopen) then
          write(iuredh,'("addwt   = T")')
          comment = 'nod pairs weighted by signal'
          if(dofits) call fithlog('ADDWT   ',.true.,comment,iufith)
        endif

        addtime = beamtime*sumwtsq/wtmax**2
        accumtim = (nz-nzero)*beamtime
        if(nodon) then
          addtime = 2.*addtime
          accumtim = 2.*accumtim
        endif
        print '(" Effective on-source time:",f7.2,"/",f7.2)', &
          addtime,accumtim
        write(iupl,'(" addtime: ",f8.2," / ",f8.2)') addtime,accumtim
        if(hdopen) then
          write(iuredh,'("accumtim= ",f8.2)') accumtim
          comment = 'accumulated on-source time per beam (s)'
          if(dofits) call fithreal('ACCUMTIM',accumtim,comment,iufith)
          if(dosum) then
            write(iusumh,'("accumtim= ",f8.2)') accumtim
            if(dofits) call fithreal('ACCUMTIM',accumtim,comment,iufish)
          endif
          write(iuredh,'("addtime = ",f8.2)') addtime
          comment = 'effective on-source time after weighting (s)'
          if(dofits) call fithreal('ADDTIME ',addtime,comment,iufith)
          if(dosum) then
            write(iusumh,'("addtime = ",f8.2)') addtime
            if(dofits) call fithreal('ADDTIME ',addtime,comment,iufish)
          endif
        endif

        scale = -1.
        if(ask) scale = 0.
        nz1 = nz+1
        print '(" Plotting collapsed shifted diffs")'
        call perspec(nx,nz1,plot,scale,iwin2)
        call waitasec(pause,.false.,plot,mx,mx,mx)

        return
      end


      subroutine adddiffs(arr,sumarr,std,nz,ierr)
!	average diffs with equal or predetermined weighting

        use ius
        use dims
        use modes

        real arr(mx,my,mp),sumarr(mx,my),std(mx,my)
        real nodpa,lores,kmirror,krot
        logical good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical hdopen
        character(60) comment

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl

        if(modeobs==MFLAT) then
          do iy = 1,ny
            do ix = 1,nx
              sumarr(ix,iy) = arr(ix,iy,1)-arr(ix,iy,2)
              std(ix,iy) = arr(ix,iy,3)
            enddo   
          enddo   
          ierr = 0
          return
        endif

        nzero = 0
        sumwt = 0.
        sumwtsq = 0.
        do iz = 1,nz
          if(wt(iz)==0.) nzero = nzero+1
          sumwt = sumwt+abs(wt(iz))
          sumwtsq = sumwtsq+wt(iz)**2
        enddo   
        if(sumwt==0.) then
          print '(" All weights = 0.  Can''t add pairs")'
          write(iupl,'(" All weights = 0.  Can''t add pairs")')
          ierr = 1
          return
        else
          stdrat = sqrt(sumwtsq)/sumwt
        endif

        do iy=1,ny
          do ix=1,nx
!	    if(illum(ix,iy)>=0) then
            if(illum(ix,iy)>0) then
              suma = 0.
              do iz=1,nz
                if(wt(iz)/=0.) then
                  suma = suma+wt(iz)*arr(ix,iy,iz)
                endif
              enddo   
              sumarr(ix,iy) = suma/sumwt
              std(ix,iy) = std(ix,iy)*stdrat
            else
              sumarr(ix,iy) = 0.
              std(ix,iy) = 0.
            endif
          enddo   
        enddo   

        inquire(unit=iuredh,opened=hdopen)
        if(hdopen) then
          write(iuredh,'("addwt   = F")')
          comment = 'nod pairs added with equal weight'
          if(dofits) call fithlog('ADDWT   ',.false.,comment,iufith)
        endif

        addtime = (nz-nzero)*beamtime
        if(nodon) addtime = 2.*addtime
        if(hdopen) then
          write(iuredh,'("addtime = ",f11.6)') addtime
          comment = 'accumulated on-source integration time'
          if(dofits) call fithreal('ACCUMTIM',addtime,comment,iufith)
          if(dosum) then
            write(iusumh,'("addtime = ",f8.2)') addtime
            if(dofits) call fithreal('ACCUMTIM',addtime,comment,iufish)
          endif
        endif

        ierr = 0
        return
      end


      subroutine extract(arr,flat,std,spec,outfile,ierr)
!	extract 1-d spectrum from 256^2 array, with optional weighting
!	along slit.  also calculate wavenumber scale and plot.

        use ius
        use dims
        use modes
        use consts
        use paths

        real arr(mx,my),flat(mx,my),std(mx,my)
        real spec(ms,mt)
        real templ(mx,my),extwt(my)
        real chisq(20),ddids(4)
        real nodpa,lores,kmirror,krot,lskrot
        real black(ms),plot(ms,8),plsm(mx),bslo(ms),bshi(ms),bsmid(ms)
        character(32) outfile,nomon
        character(40) atmofile
        character(60) comment
        character(8) yn
        logical baddata,doplot
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical good,hdopen
        logical submean,flatnorm,doextwt,sumtempl,scan,justwno
        character(32) rawfile,cardfile,lastcard,redfile,sumfile,posfile
        character(48) rawfits,cardfits,redfits,sumfits

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /tort0/ hrfl0,xdfl0,xdmrdgr,xdlrdgr,xdmr0,xdlr0,fred0, &
                xdkrot,lskrot
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /iwin/ iwin1,iwin2
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime
        common /filenames/rawfile,cardfile,lastcard,redfile,sumfile,posfile, &
                rawfits,cardfits,redfits,sumfits,ired

        save wnomin,wnomax,resolv

        if(dofits) then
          inquire(unit=iufith,opened=hdopen)
          if(.not.hdopen) open(unit=iufith,file='fits.hd',position='append')
        endif
!	atmodir should be set in mods.h
!        atmodir = '/home/texes/pipe/'	! IRTF
        lenatmo = index(atmodir,' ')-1
        doplot = (iwin1>0)
        nomon = ' '
        justwno = (ierr==-1)
        ierr = 0
!	Warning: slitoff must also be set in storescan
!	But what is it for?  It seems to make hrr that sets the dispersion
!	differ from that derived from the position of the calibration line.
        slitoff = -0.3
        if(wno0==0.) then
          wnoc = waveno0
        else
          wnoc = abs(wno0)
        endif

!	make template from array collapsed and replicated spectrally

        submean = (dosubsky &
          .and.((modeobs==MNOD).or.(modeobs==MMAP)))
        flatnorm = .false.
        sumtempl = .false.
        scan = ((modeobs==MSCAN).or.(modeobs==MMAP))
        doextwt = (.not.justwno).and.(.not.scan).and. &
          ((modeext==MCONWT).or.(modeext==MCONINT))
        if(doextwt) &
          call maketempl(arr,flat,std,templ,1, &
            submean,flatnorm,sumtempl,ierr)

      if((modeinst==MHIMED).or.(modeinst==MHILOW)) then

!	cross-dispersed modes; spectrum runs along y

        if(norder*ny>ms) then
          print '("*** Clipping off orders to fit in spec array",2i4)', &
            norder,ms/ny
          norder = ms/ny
          ns = norder*ny
        endif
        ispac = int(spacing)
        if(ispac/=nt) print '("***ispac, nt:",2i4)', ispac,nt
        if((nt+4)>mt) then
          nt = mt-4
          print '(" nt+4 > mt  Setting nt =",i4)', nt
        endif
        nt4 = nt+4
        ny2 = ny/2
        if((modeext>=MINTWT) &
          .and.((intext(1)<1).or.(intext(1)>nt) &
          .or.(intext(2)<intext(1)).or.(intext(2)>nt) &
          .or.(intext(3)<0).or.(intext(3)>nt) &
          .or.(intext(4)<intext(3)).or.(intext(4)>nt))) &
          print '(" ** Warning: extint set outside of slit range",6i4)', &
            (intext(i),i=1,4),1,nt
        if((modeext==MNODWT) &
          .and.((intext(2)<intext(1)).or.(intext(2)>nt))) &
          print '(" ** Warning: extint set outside of slit range",6i4)', &
            (intext(i),i=1,4),1,nt

!	initialize extwt over extraction interval

        do it = 1,nt
          extwt(it) = 0.
        enddo   
        if((intill(1)>0).and.(intill(2)>intill(1)) &
          .and.(intill(2)<=nt)) then
          it1 = intill(1)
          it2 = intill(2)
          it4 = intill(2)
        else
          it1 = 0
          it4 = 0
          do it = 1,nt
            if((it1==0).and.(illx(it)>0)) it1 = it
            if(illx(it)>0) it4 = it
          enddo
          it2 = it4
        endif
        if((intext(1)>0).and.(intext(2)>intext(1))) then
          if((intext(3)>0).and.(intext(4)>intext(3))) then
            jt1 = min(intext(1),intext(3))
            jt4 = max(intext(2),intext(4))
          else
            jt1 = intext(1)
            jt4 = max(intext(2),intext(4))
          endif
        else
          ntill = it4-it1+1
          jt1 = it1+ntill/8
          jt4 = it4-ntill/8
        endif
        jdt = max(1,(jt4-jt1-1)/2)
        jt2 = jt1+jdt
        jt3 = jt4-jdt
        jt5 = jt1+jdt/2
        jt6 = jt4-jdt/2
        if((jt2<jt1).or.(jt6<jt5)) then
          print '("*** Bad flat extraction interval: ",6i4)', &
            jt1,jt2,jt3,jt4,jt5,jt6
          stop
        endif

        if((modeext==MUNWT).or.(modeext==MCONWT) &
          .or.(modeext==MNONE)) then
          do it = it1,it2
            extwt(it) = 1.
          enddo   
        elseif(modeext==MNODWT) then
          if(intext(2)>intext(1)) then
            it1 = max(it1,intext(1))
            it4 = min(it4,intext(2))
          elseif(intext(2)<intext(1)) then
            it1 = max(it1,intext(2))
            it4 = min(it4,intext(1))
          endif
          idt = (it4-it1-1)/2
          it2 = it1+idt
          it3 = it4-idt
          if((date>6.0).neqv.(intext(2)<intext(1))) then
            do it = it1,it2
              extwt(it) = -1.
            enddo   
            do it = it3,it4
              extwt(it) = 1.
            enddo   
          else
            do it = it1,it2
              extwt(it) = 1.
            enddo   
            do it = it3,it4
              extwt(it) = -1.
            enddo   
          endif
        elseif((modeext==MINTWT).or.(modeext==MAVGWT)) then
          if(intext(1)*intext(2)==0) then
            print '("*** extint not set." &
              "  Can''t extract with extmode = int")'
            ierr = 9
            return
          elseif(intext(3)*intext(4)==0) then
            it1 = max(it1,intext(1))
            it2 = intext(2)
            if((it2>it4).or.(it2<it1)) it2 = it4
            it3 = 0
            it4 = 0
            do it = it1,it2
              extwt(it) = 1.
            enddo   
          else
            it1 = max(it1,intext(1))
            it2 = intext(2)
            if((it2>it4).or.(it2<it1)) it2 = it4
            it3 = intext(3)
            it4 = min(it4,intext(4))
            if((it3>0).and.(it4>it3)) then
              it4 = min(it4,nt)
              do it = it3,it4
                extwt(it) = -1.
              enddo   
            endif
            do it = it1,it2
              extwt(it) = 1.
            enddo   
          endif
        elseif(modeext==MCONINT) then
          if(intext(1)*intext(2)==0) then
            print '("*** extint not set." &
              "  Can''t extract with extmode = conint")'
            ierr = 9
            return
          endif
          it1 = max(it1,intext(1))
          it2 = intext(2)
          if((it2>it4).or.(it2<it1)) it2 = it4
          it3 = intext(3)
          it4 = min(it4,intext(4))
          if((it3>0).and.(it4>it3)) then
            it4 = min(it4,nt)
            do it = it3,it4
              extwt(it) = 1.
            enddo   
          endif
          do it = it1,it2
            extwt(it) = 1.
          enddo   
        endif
        if((it3>0).and.(it3<it1)) it1 = it3
        it4 = max(it2,it4)

!	multiply templ by extwt and normalize to preserve flux

!	convert to Jy
        if(modeext==MAVGWT) then
          omega = it4-it1+1
        else
          omega = omegap/2.998e-13
          if(nodon) omega = omega/2.
        endif
        sumsum = 0.
        avgwt = 0.
        do iorder = 1,norder
          ix0 = (iorder-1)*ispac
          if(doextwt) then
            iymx = ny/2
            sumtmx = 0.
            do iy = 1,ny
              sumt = 0.
              do it = it1,it4
                if(extwt(it)/=0.) then
                  ix = ix0+it
                  if((ix<1).or.(ix>nx)) cycle
                  templa = abs(templ(ix,iy)*extwt(it))
                  sumt = sumt+templa
                endif
              enddo   
              if(sumt>sumtmx) then
                iymx = iy
                sumtmx = sumt
              endif
            enddo
            isum = 0
            sumt = 0.
            if(dosubsky) then
              do it = it1,it4
                if(extwt(it)/=0.) then
                  ix = ix0+it
                  if((ix<1).or.(ix>nx)) cycle
                  templw = templ(ix,iymx)*extwt(it)
                  isum = isum+1
                  sumt = sumt+templw
                endif
              enddo   
              avgwt = sumt/isum
            else
              avgwt = 0.
            endif
            sumt = 0.
            sumt2 = 0.
            do it = it1,it4
              if(extwt(it)/=0.) then
                ix = ix0+it
                if((ix<1).or.(ix>nx)) cycle
                templw = templ(ix,iymx)*extwt(it)-avgwt
                sumt = sumt+abs(templw)
                sumt2 = sumt2+templw**2
              endif
            enddo   
            wtnorm = omega*sumt/sumt2
            if(0.*wtnorm/=0.) then
              print '("*** In extract wtnorm =",e10.3)', wtnorm
              wtnorm = 1.
            endif
            do it = 1,nt
              extwtit = extwt(it)
              ix = ix0+it
              if((ix<1).or.(ix>nx)) cycle
              do iy = 1,ny
                templ(ix,iy) = (templ(ix,iy)*extwtit-avgwt)*wtnorm
              enddo   
            enddo   
            sumsum = sumsum+sumt2*wtnorm**2
          else
            wtnorm = 1.
            do it = it1,it4
              ix = ix0+it
              if((ix<1).or.(ix>nx)) cycle
              if(illx(ix)>=1) then
                templi = omega*extwt(it)
              else
                templi = 0.
              endif
              sumt2 = sumt2+templi**2
              do iy = 1,ny
                templ(ix,iy) = templi
              enddo   
            enddo   
            sumsum = sumsum+sumt2
          endif
        enddo   
        lx = ix+1
        do ix = lx,nx
          do iy = 1,ny
            templ(ix,iy) = 0.
          enddo   
        enddo   
        sumt2 = sumsum/norder
        if(sumt2==0.) sumt2 = 1.

!	weight and sum spatially within orders

        illsum = 0
        do iorder = 1,norder
          ix0 = (iorder-1)*ispac
          do it = it1,it4
            ix = ix0+it
            if((ix>=1).and.(ix<=nx).and.(illx(ix)==1)) &
              illsum = illsum+1
          enddo   
        enddo   
!	max number of non-illuminated pixels in an order
!        illmax = (it4-it1+1)-(illsum/norder)+2
        illmax = min(2,illsum/(8*norder))
        specmax = 1.

!	make wno scale, noise and atmo spectra

!	dlnw = pixelwd/(2.*abs(hrr)*hrfl)
  280   dlnw = pixelwd/(2.*abs(hrr)*(1.-slitoff/20.)*hrfl)
        nsumsum = 0
        do iorder = 1,norder
          ix0 = (iorder-1)*ispac
          is0 = ny*(iorder-1)
          dw = 0.5/(sqrt(hrr**2/(1.+hrr**2))*hrdgr)
          wnoi = wnoc+(iorder-(norder+1)/2)*dw
          do iy = 1,ny
            y = iy-ny2
            is = is0+iy
            sum = 0.
            sumvar = 0.
            sumcard = 0.
            sumblk = 0.
            sumsky = 0.
            sumwt = 0.
            sumwt2 = 0.
            sumt2 = 0.
            sumlo = 0.
            sumhi = 0.
            summid = 0.
            nsumlo = 0
            nsumhi = 0
            nsummid = 0
            illsum = 0
            do it = it1,it4
              jt = it+4
              ix = ix0+it
              if((ix>=1).and.(ix<=nx).and.(illx(ix)==1)) then
                if(illum(ix,iy)==1) then
                  templi = templ(ix,iy)
                  templa = abs(templi)
                  if(justwno) then
                    spec(is,jt) = cards(ix,iy,2)
                  else
                    spec(is,jt) = arr(ix,iy)
                    sum = sum+templi*arr(ix,iy)
                    if(scan) then
                      sumvar = sumvar+templa*std(ix,iy)
                    else
                      sumvar = sumvar+(templi*std(ix,iy))**2
                    endif
                  endif
                  if((modeobs==MSTARE).and.(modecard==MNONE)) then
                    sumcard = sumcard+templa*arr(ix,iy)
                    sumblk = sumblk+templa
                  elseif(cards(ix,iy,2)/=0.) then
                    sumcard = sumcard+templa*cards(ix,iy,2)
                    sumblk = sumblk+templa*cards(ix,iy,1)
                    sumsky = sumsky+templa*cards(ix,iy,3)
                    if((it>=jt1).and.(it<=jt2)) then
                      sumlo = sumlo+cards(ix,iy,2)
                      nsumlo = nsumlo+1
                    elseif((it>=jt3).and.(it<=jt4)) then
                      sumhi = sumhi+cards(ix,iy,2)
                      nsumhi = nsumhi+1
                    endif
                    if((it>=jt5).and.(it<=jt6)) then
                      summid = summid+cards(ix,iy,2)
                      nsummid = nsummid+1
                      nsumsum = nsumsum+1
                    endif
                  endif
                  sumwt = sumwt+templa
                  sumwt2 = sumwt2+templa**2
!                elseif((illum(ix,iy)==-1).or.(illsum>illmax)) then
                elseif(illsum>illmax) then
                  do jt = 5,nt4
                    spec(is,jt) = 0.
                  enddo   
                  sum = 0.
                  sumvar = 0.
                  sumcard = 0.
                  sumblk = 0.
                  sumwt2 = 0.
                  sumlo = 0.
                  sumhi = 0.
                  summid = 0.
                  exit
                else
                  illsum = illsum+1
                  spec(is,jt) = 0.
                endif
                sumt2 = sumt2+templa**2
              else
                spec(is,jt) = 0.
              endif
            enddo   
            do jt = nt4+1,mt
              spec(is,jt) = 0.
            enddo   
            renorm = 1.
            if((doextwt).and.(sumwt2>0.)) renorm = sumt2/sumwt2
            if(sumwt<=0.) sumwt = 1.
            spec(is,1) = wnoi*sexp(dlnw*y)
            spec(is,2) = sum*renorm
            if(scan) then
              spec(is,3) = sumvar/sumwt
            else
              spec(is,3) = sqrt(sumvar)*renorm
            endif
            spec(is,4) = sumcard/sumwt
!            if((modecard==MSKY).or.(modecard==MOBJ)) then
!              black(is) = sumsky/sumwt
!            else
              black(is) = sumblk/sumwt
!            endif
            specmax = max(specmax,sum*renorm)
            if((nsumlo>jdt/2).and.(nsumhi>jdt/2)) then
              bslo(is) = sumlo/nsumlo
              bshi(is) = sumhi/nsumhi
            else
              bslo(is) = 0.
              bshi(is) = 0.
            endif
            if(nsummid>jdt/2) then
              bsmid(is) = summid/nsummid
            else
              bsmid(is) = 0.
            endif
          enddo   
        enddo   
        if((modecard/=MNONE).and.(nsumsum==0)) then
          print '("*** All cards = 0 in extract?")'
          stop
        endif

        ns = ny*norder
        do it = 1,mt
          do is = ns+1,ms
            spec(is,it) = 0.
          enddo   
        enddo   

        smax = 0.
        smin = 0.
        nplt = norder
        nsmoo = 2

        do iorder = 1,norder
          is0 = (iorder-1)*ny
          do iy = 1,ny
            is = is0+iy
            plot(is,1) = spec(is,1)
            if(modecard==MNONE) then
              plot(is,2) = spec(is,2)
            else
              plot(is,2) = bsmid(is)
            endif
            plot(is,3) = 0.
            plot(is,4) = 0.
            plot(is,5) = 0.
            plot(is,7) = bslo(is)
            plot(is,8) = bshi(is)
            if(plot(is,2)>smax) smax = plot(is,2)
          enddo   
        enddo   
        smin = 0.5*min(smax,1.)
        if(ask) then
          scale = 0.
        elseif(modecard==MNONE) then
          scale = smax
        else
          scale = 1.
        endif
        if(wno0==0.) then
          print '(" Plotting reference spectrum")'
          call plots(ms,ny,nplt,plot,scale,nomon,iwin2)
  270     print '(" Enter reference wavenumber and order ",$)'
          read *, wnoref,iorder
          if((wnoref<=0.).or.(iorder<1)) goto 299
          is0 = ny*(iorder-1)
          do iy = 1,ny
            is = is0+iy
            plot(iy,1) = spec(is,1)
            if(modecard==MNONE) then
              plot(iy,2) = spec(is,2)
            else
              plot(iy,2) = bsmid(is)
            endif
          enddo   
          call plots(ms,ny,1,plot,scale,nomon,iwin2)
          print '(" Click cursor on reference line")'
          call getcurs(x,y,ikey,ierr)
!	 w = x
          if((modecard==MSKY).or.(modecard==MNONE)) then
            w = plmax(plot,ms,x,nsmoo)
          else
            w = plmin(plot,ms,x,nsmoo)
          endif
          dw = 0.5/(sqrt(hrr**2/(1.+hrr**2))*hrdgr)
          wnoi = wnoc+(iorder-(norder+1)/2)*dw
          wnoc = wnoc+(wnoref-w)*wnoi/w
          wno0 = wnoc
          print '(" wno0 = ",f10.3)', wno0
          write(iupl,'("wno0 = ",f10.3)') wno0
          if(modeinst==MHIMED) then
            xdang = echelle
          else
            xdang = lores
          endif
!	this might not be right
          xdorder = nint(2.0*xddgr*wno0*sin(xdang/DEGRAD))
!	and wno0 isn''t really at center of array
          xdang0 = DEGRAD*asin(xdorder/(2.0*xddgr*wno0))
          if(modeinst==MHIMED) then
            xdmr0 = echelle+xdmr0-xdang0
            echelle = xdang0
            print '(" mr0 = ",f10.3)', xdmr0
            write(iupl,'(" mr0 = ",f10.3)') xdmr0
          else
            xdlr1 = xdlr0
            xdlr0 = lores+xdlr0-xdang0
            lores = xdang0
            print '(" lr0 = ",f10.3)', xdlr0
            write(iupl,'(" lr0 = ",f10.3)') xdlr0
          endif
          hrorder = 2.0*hrdgr*wno0*sin(atan(abs(hrr)))
          hrr0 = tan(asin(nint(hrorder)/(2.0*hrdgr*wno0)))
          print '(" Calculated echelon order, R# = ",f8.2,f8.3)', &
            hrorder,hrr0
          hrrb = 10.0-200./wno0
          if((hrr>0.).and.(abs(hrr0-hrrb)>400./wno0)) then
            hrorder = nint(2.0*hrdgr*wno0*sin(atan(hrrb)))
            hrr0 = tan(asin(hrorder/(2.0*hrdgr*wno0)))
            print '(" Doesn''t sound right.  Instead let''s try ", &
              f8.2,f8.3)', hrorder,hrr0
          endif
          write(iupl,'(" Calculated echelon order, R# = ",f8.2,f8.3)') &
            hrorder,hrr0
          print '(" Current R# = ",f8.3, &
            "  Change hrr or smoothing? (n) ",$)', hrr
          read '(a8)', yn
          if(index(yn,'y')>0) then
            print '(" Enter hrr, nsmooth ",$)'
            read *, hrr,nsmoo
!	    dlnw = pixelwd/(2.*abs(hrr)*hrfl)
            dlnw = pixelwd/(2.*abs(hrr)*(1.-slitoff/20.)*hrfl)
            dw = 0.5/(sqrt(hrr**2/(1.+hrr**2))*hrdgr)
          endif
          do iorder = 1,norder
            is0 = (iorder-1)*ny
            wnoi = wnoc+(iorder-(norder+1)/2)*dw
            do iy = 1,ny
              y = iy-ny2
              is = is0+iy
              spec(is,1) = wnoi*sexp(dlnw*y)
              plot(is,1) = spec(is,1)
!              plot(is,2) = spec(is,4)
              plot(is,2) = bsmid(is)
            enddo   
          enddo   
          call plots(ms,ny,nplt,plot,scale,nomon,iwin2)
          if(index(yn,'y')>0) goto 270
          if(scan) ierr = 1
        elseif((wno0<0.).and.(modecard/=MNONE)) then
!	read atmo file and interpolate onto spec wno grid
          wno0 = abs(wno0)
          if(wno0<500) then
            atmofile = 'nocando'
          elseif(wno0<700) then
            atmofile = atmodir(1:lenatmo)//'Q.mod'
          elseif(wno0<1000) then
            atmofile = atmodir(1:lenatmo)//'NL.mod'
          elseif(wno0<1400.) then
            atmofile = atmodir(1:lenatmo)//'NS.mod'
          elseif(wno0<1800.) then
            atmofile = 'nocando'
          elseif(wno0<2300.) then
            atmofile = atmodir(1:lenatmo)//'M.mod'
          else
            atmofile = 'nocando'
          endif
          open(unit=iuatmo,file=atmofile,access='sequential',iostat=jerr)
          if(jerr/=0) then
            print '(" Error opening ",a32)', atmofile
            print '(" May need to change atmodir in mods.h")'
            print '(" Hit RETURN to continue or ^C to quit")'
            read '(a8)', yn
            goto 299
          endif
          dwnoat = 0.001
          airpow = max(1.,min(3.,airmass))
          do is = 1,ns
            wnois = spec(is,1)
!	if needed, step back to last atmo wno before wno(is)
            do
              read(unit=iuatmo,fmt=*,iostat=jerr) &
                wnoa,trhrlw,trhrhw,bshrlw,bshrhw,bslrlw,bslrhw
              if(jerr/=0) then
                print '(" Error reading atmo model",2f8.3)', wnois,wnoa
                goto 299
              endif
              if(wnoa<wnois-dwnoat) exit
              backspace(unit=iuatmo)
              backspace(unit=iuatmo)
            enddo
!	step to last atmo wno before wno(is)
            do
              read(unit=iuatmo,fmt=*,iostat=jerr) &
                wnoa,trhrlw,trhrhw,bshrlw,bshrhw,bslrlw,bslrhw
              if(jerr/=0) then
                print '(" Error reading atmo model",2f8.3)', wnois,wnoa
                goto 299
              endif
              if(wnoa>wnois-dwnoat) exit
            enddo
            bshia = max(0.,bshrhw)
            bsloa = max(0.,bshrlw)
            trhia = max(0.,trhrhw)
            trloa = max(0.,trhrlw)
!	read atmo wno after wno(is)
            read(unit=iuatmo,fmt=*,iostat=jerr) &
                wnoa,trhrlw,trhrhw,bshrlw,bshrhw,bslrlw,bslrhw
            bshib = max(0.,bshrhw)
            bslob = max(0.,bshrlw)
            trhib = max(0.,trhrhw)
            trlob = max(0.,trhrlw)
            diat = (wnois-wnoa)/dwnoat
            bshis = max(0.,bshia*(1.0-diat)+bshib*diat)
            bslos = max(0.,bsloa*(1.0-diat)+bslob*diat)
            trhis = max(0.,trhia*(1.0-diat)+trhib*diat)
            trlos = max(0.,trloa*(1.0-diat)+trlob*diat)
            plot(is,3) = bshis**airpow
            if((bshis<0.).and.(bslos>0.)) then
              plot(is,4) = bslos**(1.75*airpow)
            else
              plot(is,4) = (bshis*bslos)**(airpow/2.)
            endif
            plot(is,5) = bslos**airpow
!	this may be overwritten by sky emission
            plot(is,6) = (trhis/2.+trlos/2.)**airpow
            plot(is,7) = bslo(is)
            plot(is,8) = bshi(is)
          enddo
          close(unit=iuatmo)
          print '(" Plotting blk-sky, models")'
          call plot4(ms,ny,nplt,plot,scale,nomon,iwin2)
          if(testrun) then
            call waitasec(.true.,.false.,plot,ms,ms,ny)
          else
            call waitasec(pause,.false.,plot,ms,ms,ny)
          endif
!	subtract smoothed spectra
          do iorder = 1,norder
            do ipl = 2,5
              do iy = 1,ny
                is = ny*(iorder-1)+iy
                ky = min(32,iy-1,ny-iy)
                yk = ky+1.
                sumpl = 0.
                sumwt = 0.
                do jdy = -ky,ky
                  js = is+jdy
                  if(plot(js,2)/=0.) then
                    wtj = 1.-abs(jdy)/yk
                    sumpl = sumpl+wtj*plot(js,ipl)
                    sumwt = sumwt+wtj
                  endif
                enddo
                if(sumwt>0.) then
                  plsm(iy) = sumpl/sumwt
                else
                  plsm(iy) = 0.
                endif
              enddo
              do iy = 1,ny
                is = ny*(iorder-1)+iy
                if(plot(is,2)/=0.) then
                  plot(is,ipl) = plot(is,ipl)-plsm(iy)
                else
                  plot(is,ipl) = 0.
                endif
              enddo
            enddo
          enddo
!	find shift between blk-sky and atmo model
          if(norder>4) then
            nds = 8
            norder1 = 2
            norder2 = norder-1
          else
            nds = 7
            norder1 = 1
            norder2 = norder
          endif
          idsmin = -nds
          chisqmin = 1.0e+06
          do ids = -nds,nds
            if(abs(ids)==8) then
              jds = sign(256,ids)
            else
              jds = ids
            endif
            sumerr3 = 0.
            sumsq3 = 0.
            sumerr4 = 0.
            sumsq4 = 0.
            sumerr5 = 0.
            sumsq5 = 0.
            sumn = 0.
            do iorder = norder1,norder2
              do iy = 8,ny-7
                is = ny*(iorder-1)+iy
                if(plot(is,2)*plot(is-1,2)*plot(is+1,2)>0.) then
                  dbs = plot(is,2)-plot(is+jds,3)
                  sumerr3 = sumerr3+dbs
                  sumsq3 = sumsq3+dbs**2
                  dbs = plot(is,2)-plot(is+jds,4)
                  sumerr4 = sumerr4+dbs
                  sumsq4 = sumsq4+dbs**2
                  dbs = plot(is,2)-plot(is+jds,5)
                  sumerr5 = sumerr5+dbs
                  sumsq5 = sumsq5+dbs**2
                  sumn = sumn+1.
                endif
              enddo
            enddo
            if(sumn==0.) sumn = 1.
            chisq3 = sumsq3-sumerr3**2/sumn
            chisq4 = sumsq4-sumerr4**2/sumn
            chisq5 = sumsq5-sumerr5**2/sumn
            slope = (chisq5-chisq3)/2.
            curve = max(0.,(chisq3-2.*chisq4+chisq5))
            if(slope==0.) then
              print '(" Bad chisq for pixel shift ",i2,5f9.3)', &
                ids,chisq3,chisq4,chisq5,slope,curve
              chisq(ids+9) = 100.
              cycle
            endif
            dmin = -slope/max(abs(slope),curve)
            chisqi = chisq4+slope*dmin+curve*dmin**2/2.
            chisq(ids+9) = chisqi
            if(chisqi<chisqmin) then
              idsmin = ids
              chisqmin = chisqi
            endif
          enddo
          if(abs(idsmin)==8) then
            ishift = idsmin/8
            print '("*** Best wno order shift = ",i2)', ishift
            print '(" Shift by ",i2," order? (y) "$)', ishift
            read '(a8)', yn
            if(index(yn,'n')==0) then
              wno0 = wno0+0.6623*ishift
              wnoc = wno0
              wno0 = -wno0
              goto 280
            endif
          endif
          if(abs(idsmin)>2) then
            print '("*** Best wno pixel shift = ",i2)', idsmin
  285       print '(" Shift by p(ixel), w(no), o(rder), n(one),", &
              " or s(kip), k(bd) ",$)'
            read '(a8)', yn
            if(index(yn,'n')>0) then
              dids = 0.
              goto 290
            elseif(index(yn,'s')>0) then
              ierr = 8
              return
            elseif(index(yn,'k')>0) then
              ierr = 9
              return
            else
              print '(" Enter shift ", &
                "(positive to shift observed spectrum blue) ",$)'
              read(*,*,iostat=ierr) shift
              if(ierr>0) goto 285
              if(shift==0.) then
                dids = 0.
                goto 290
              elseif(index(yn,'p')>0) then
                wno0 = wno0+shift*(plot(192,1)-plot(64,1))/128.
              elseif(index(yn,'w')>0) then
                wno0 = wno0+shift
              elseif(index(yn,'o')>0) then
                wno0 = wno0+0.6623*shift
              else
                goto 285
              endif
              wnoc = wno0
              wno0 = -wno0
              goto 280
            endif
          else
            chisq1 = chisq(idsmin+8)
            chisq2 = chisq(idsmin+9)
            chisq3 = chisq(idsmin+10)
            dids = ((chisq1-chisq3)/2.)/(chisq1+chisq3-2.*chisq2)
            if((abs(idsmin)>2).or.(abs(dids)>1.)) &
              print '(8f10.3/7f10.3,i6)', (chisq(ids),ids=2,16),idsmin
            dids = idsmin+min(1.,max(-1.,dids))
          endif
  290     wno0 = abs(wno0)+dids*(plot(192,1)-plot(64,1))/128.
          hrorder = 2.0*hrdgr*wno0*sin(atan(abs(hrr)))
          hrr0 = tan(asin(nint(hrorder)/(2.0*hrdgr*wno0)))
          print '(" Pixel shift, new wno0, hrr:",3f10.4)', dids,wno0,hrr0
          write(iupl,'(" Pixel shift, new wno0, hrr:",3f10.4)') &
            dids,wno0,hrr0
          if(hrr>0.) then
            hrrb = 10.0-200./wno0
            if(abs(hrr0-hrrb)>400./wno0) then
              hrorderb = nint(2.0*hrdgr*wno0*sin(atan(hrrb)))
              hrrb = tan(asin(hrorderb/(2.0*hrdgr*wno0)))
              print '(" hrr doesn''t sound right.  Instead try ", &
                f8.2,f8.3,"? "$)', hrorderb,hrrb
              read '(a8)', yn
              if(index(yn,'n')==0) then
                hrorder = hrorderb
                hrr0 = hrrb
              endif
            endif
            write(iupl,'(" Calculated echelon order, R# = ",f8.2,f8.3)') &
              hrorder,hrr0
            hrr = hrr0
          else
            print '(" hrr is fixed at ",f8.3)', -hrr
          endif
          dlnw = pixelwd/(2.*abs(hrr)*(1.-slitoff/20.)*hrfl)
          dw = 0.5/(sqrt(hrr**2/(1.+hrr**2))*hrdgr)
          do iorder = 1,norder
            is0 = ny*(iorder-1)
            wnoi = wno0+(iorder-(norder+1)/2)*dw
            do iy = 1,ny
              y = iy-ny2
              is = is0+iy
              spec(is,1) = wnoi*sexp(dlnw*y)
            enddo   
          enddo   
!	check slitrot
          do ipl = 7,8
            do iorder = 1,norder
              do iy = 1,ny
                is = ny*(iorder-1)+iy
                ky = min(16,iy-1,ny-iy)
                yk = ky+1.
                sumpl = 0.
                sumwt = 0.
                do jdy = -ky,ky
                  js = is+jdy
                  if(plot(js,ipl)/=0.) then
                    wtj = 1.-abs(jdy)/yk
                    sumpl = sumpl+wtj*plot(js,ipl)
                    sumwt = sumwt+wtj
                  endif
                enddo
                if(sumwt>0.75*yk) then
                  plsm(iy) = sumpl/sumwt
                else
                  plsm(iy) = 0.
                endif
              enddo
              do iy = 1,ny
                is = ny*(iorder-1)+iy
                if((plsm(iy)>0.).and.(plot(is,ipl)/=0.)) then
                  plot(is,ipl) = plot(is,ipl)-plsm(iy)
                else
                  plot(is,ipl) = 0.
                endif
              enddo
            enddo
            do ids = idsmin-1,idsmin+1
              sumerr3 = 0.
              sumsq3 = 0.
              sumerr4 = 0.
              sumsq4 = 0.
              sumerr5 = 0.
              sumsq5 = 0.
              sumn = 0.
              do iorder = 1,norder
                do iy = 8,ny-7
                  is = ny*(iorder-1)+iy
                  if(plot(is,ipl)*plot(is-1,ipl)*plot(is+1,ipl)/=0.) then
                    dbs = plot(is,ipl)-plot(is+ids,3)
                    sumerr3 = sumerr3+dbs
                    sumsq3 = sumsq3+dbs**2
                    dbs = plot(is,ipl)-plot(is+ids,4)
                    sumerr4 = sumerr4+dbs
                    sumsq4 = sumsq4+dbs**2
                    dbs = plot(is,ipl)-plot(is+ids,5)
                    sumerr5 = sumerr5+dbs
                    sumsq5 = sumsq5+dbs**2
                    sumn = sumn+1.
                  endif
                enddo
              enddo
              if(sumn==0.) sumn = 1.
              chisq3 = sumsq3-sumerr3**2/sumn
              chisq4 = sumsq4-sumerr4**2/sumn
              chisq5 = sumsq5-sumerr5**2/sumn
              slope = (chisq5-chisq3)/2.
              curve = max(0.,(chisq3-2.*chisq4+chisq5))
              dmin = -slope/max(abs(slope/2.),curve)
              chisqi = chisq4+slope*dmin+curve*dmin**2/2.
              chisq(ids-idsmin+2) = chisqi
            enddo
            chisq1 = chisq(1)
            chisq2 = chisq(2)
            chisq3 = chisq(3)
            ihilo = 2*ipl-13
            ddids(ihilo) = ((chisq1-chisq3)/2.)/(chisq1+chisq3-2.*chisq2)
          enddo
          ddids(2) = dids
          if(abs(ddids(1)-ddids(3))>0.4) then
            print '("*** Spectral shift along slit: ",3f10.3)', &
              (ddids(ihilo),ihilo=1,3)
          else
            print '(" Spectral shift along slit: ",3f10.3)', &
              (ddids(ihilo),ihilo=1,3)
          endif
!	check xd grating0
          if(modeinst==MHIMED) then
            xdang = echelle
          else
            xdang = lores
          endif
          xdorder = nint(2.0*xddgr*wno0*sin(xdang/DEGRAD))
          xdang0 = DEGRAD*asin(xdorder/(2.0*xddgr*wno0))
          if(modeinst==MHIMED) then
            xdmr0 = echelle+xdmr0-xdang0
            echelle = xdang0
            print '(" mr0 = ",f10.3)', xdmr0
            write(iupl,'(" mr0 = ",f10.3)') xdmr0
          else
            xdlr0 = lores+xdlr0-xdang0
            lores = xdang0
            print '(" lr0 = ",f10.3)', xdlr0
            write(iupl,'(" lr0 = ",f10.3)') xdlr0
          endif
          if(rdfits) then
            write(iuwno,'(a32,2f8.3)') rawfits,wno0,hrr
          else
            write(iuwno,'(a24,2f8.3)') rawfile,wno0,hrr
          endif
        else
          print '(" Plotting reference spectrum")'
          call plot4(ms,ny,nplt,plot,scale,nomon,iwin2)
          call waitasec(pause,.false.,plot,ms,ms,ny)
        endif

  299   if(justwno.or.(modeext==MNONE)) goto 289
        sumsn = 0.
        sumn= 0.
        smax = -1.
        dopp = (1.+radvel/CLIGHT)
        do iorder = 1,norder
          do ip = 1,ny
            is = (iorder-1)*ny+ip
            plot(is,1) = dopp*spec(is,1)
            plot(is,2) = spec(is,2)
            plot(is,3) = spec(is,3)
            plot(is,4) = spec(is,4)
            plot(is,5) = black(is)
            plot(is,6) = black(is)*(1.0-spec(is,4))
            if(spec(is,4)>smin) then
              if(spec(is,2)>smax) smax = spec(is,2)
              if(spec(is,3)>0.) then
                sumsn = sumsn+spec(is,2)/spec(is,3)
                sumn = sumn+1.
              endif
            endif
          enddo   
        enddo   
        sumsn = sumsn/sumn
        print '(" Mean S/N in extracted spectrum: ",f10.2)', sumsn
        write(iupl,'(" Mean S/N in extracted spectrum: ",f10.2)') sumsn
        scale = -1.
        if(ask) scale = 0.
        print '(" Plotting extracted spectrum")'
        if(domon) then
          call plots(ms,ny,nplt,plot,scale,outfile,iwin2)
        else
          call plots(ms,ny,nplt,plot,scale,nomon,iwin2)
        endif
        call waitasec(pause,.false.,plot,ms,ms,ny)

!	dlnw = pixelwd/(2.*abs(hrr)*hrfl)
        dlnw = pixelwd/(2.*abs(hrr)*(1.-slitoff/20.)*hrfl)
        wnomin = spec(1,1)
        wnomax = spec(ns,1)
        dvpix = -CLIGHT*dlnw
        rslit = pltscl/(dlnw*slitwid)
        rdiff = 18.0*abs(hrr)*wno0
        resolv = 1.0/sqrt(1.0/rslit**2+1.0/rdiff**2)
        inquire(unit=iuredh,opened=hdopen)
        if(hdopen) then
          write(iuredh,'("spacing = ",f10.6)') spacing
          write(iuredh,'("xorder1 = ",f10.6)') xorder1
          write(iuredh,'("hrr     = ",f10.4)') abs(hrr)
          write(iuredh,'("wno0    = ",f10.4)') wno0
          write(iuredh,'("wnomin  = ",f10.4)') wnomin
          write(iuredh,'("wnomax  = ",f10.4)') wnomax
          write(iuredh,'("dvpix   = ",f10.4)') dvpix
          write(iuredh,'("resolv  = ",f10.3)') resolv
        endif
        inquire(unit=iufith,opened=hdopen)
        if(dofits.and.hdopen) then
          comment = 'order separation in pixels'
          call fithreal('SPACING ',spacing,comment,iufith)
          comment = 'first pixel of order 1'
          call fithreal('XORDER1 ',xorder1,comment,iufith)
          comment = 'R number of echelon grating'
          call fithreal('HRR     ',abs(hrr),comment,iufith)
          comment = 'wavenumber at pixel 128 in int(central order)'
          call fithreal('WNO0    ',wno0,comment,iufith)
          comment = 'minimum wavenumber (cm-1)'
          call fithreal('WNO_MIN ',wnomin,comment,iufith)
          comment = 'maximum wavenumber (cm-1)'
          call fithreal('WNO_MAX ',wnomax,comment,iufith)
          comment = 'velocity interval between pixels (km/s)'
          call fithreal('DVPIX   ',dvpix,comment,iufith)
          comment = 'approximate resolving power'
          call fithreal('RESOLV  ',resolv,comment,iufith)
        elseif(dofits) then
          print '("*** Can''t write WNO0, etc.  fits.hd not open.")'
        endif
        if(dosum.and.(sumtime==0.))  then
          inquire(unit=iusumh,opened=hdopen)
          if(hdopen) then
            write(iusumh,'("spacing = ",f10.6)') spacing
            write(iusumh,'("xorder1 = ",f10.6)') xorder1
            write(iusumh,'("hrr     = ",f10.4)') abs(hrr)
            write(iusumh,'("wno0    = ",f10.4)') wno0
            write(iusumh,'("wnomin  = ",f10.4)') wnomin
            write(iusumh,'("wnomax  = ",f10.4)') wnomax
            write(iusumh,'("dvpix   = ",f10.4)') dvpix
            write(iusumh,'("resolv  = ",f10.3)') resolv
          endif
          inquire(unit=iufish,opened=hdopen)
          if(dofits.and.hdopen) then
            comment = 'order separation in pixels'
            call fithreal('SPACING ',spacing,comment,iufish)
            comment = 'first pixel of order 1'
            call fithreal('XORDER1 ',xorder1,comment,iufish)
            comment = 'R number of echelon grating'
            call fithreal('HRR     ',abs(hrr),comment,iufish)
            comment = 'wavenumber at pixel 128 in int(central order)'
            call fithreal('WNO0    ',wno0,comment,iufish)
            comment = 'minimum wavenumber (cm-1)'
            call fithreal('WNO_MIN ',wnomin,comment,iufish)
            comment = 'maximum wavenumber (cm-1)'
            call fithreal('WNO_MAX ',wnomax,comment,iufish)
            comment = 'velocity interval between pixels (km/s)'
            call fithreal('DVPIX   ',dvpix,comment,iufish)
            comment = 'approximate resolving power'
            call fithreal('RESOLV  ',resolv,comment,iufish)
          elseif(dofits) then
            print '("*** Can''t write WNO0, etc.  fitsum.hd not open.")'
          endif
        endif
  289   continue

      elseif((modeinst==MMED).or.(modeinst==MLOW)) then

!	long-slit mode; spectrum runs along x
!	weight in proportion to S/N, summed spectrally

        nx2 = nx/2
        ny2 = ny/2

!	initialize extwt over extraction interval
        do iy = 1,ny
          extwt(iy) = 0.
        enddo   
 
        if(modeext==MNODWT) then
          iy1 = 0
          iy4 = 0
          do iy = 1,ny
            if(illx(iy)==1) then
              extwt(iy) = 1.
              if(iy1==0) iy1 = iy
              iy4 = iy
            else
              extwt(iy) = 0.
            endif
          enddo   
          iy2 = (iy1+iy4-1)/2
          iy3 = iy4-iy2+iy1
          if(iy3>iy2+1) extwt(iy2+1) = 0.
          do iy = iy1,iy2
            extwt(iy) = -1.
          enddo   
        elseif((modeext==MUNWT).or.(modeext==MCONWT).or.(modeext==MNONE) &
          .or.(intext(1)<=0).or.(intext(2)<intext(1))) then
          iy1 = 0
          iy4 = 0
          do iy = 1,ny
            if(illx(iy)==1) then
              extwt(iy) = 1.
              if((iy1==0).or.((iy>1).and.(illx(iy-1)==0).and.(iy<ny/2))) &
                iy1 = iy
              iy4 = iy
            else
              extwt(iy) = 0.
            endif
          enddo   
          if(iy4-iy1>80) then
            iy1 = iy1+8
            iy4 = iy4-8
          endif
          iy2 = iy4
        elseif((modeext==MINTWT).or.(modeext==MAVGWT)) then
          iy3 = intext(3)
          iy4 = min(ny,intext(4))
          if((iy3>0).and.(iy4>=iy3)) then
            do iy = iy3,iy4
              extwt(iy) = -1.
            enddo   
          endif
          iy1 = max(1,intext(1))
          iy2 = min(ny,intext(2))
          if(iy2<iy1) iy2 = ny
          do iy = iy1,iy2
            extwt(iy) = 1.
          enddo   
        elseif(modeext==MCONINT) then
          iy3 = intext(3)
          iy4 = min(ny,intext(4))
          if((iy3>0).and.(iy4>=iy3)) then
            do iy = iy3,iy4
              extwt(iy) = 1.
            enddo   
          endif
          iy1 = max(1,intext(1))
          iy2 = intext(2)
          if((iy2<iy1).or.(iy2>ny)) iy2 = ny
          do iy = iy1,iy2
            extwt(iy) = 1.
          enddo   
        endif
        iy4 = max(iy2,iy4)

!	multiply templ by extwt and normalize to preserve flux
 
!	convert to Jy
        if(modeext==MAVGWT) then
          omega = iy4-iy1+1
        else
          omega = omegap/2.998e-13
          if(nodon) omega = omega/2.
        endif
        sumsum = 0.
        if(doextwt) then
          sumt = 0.
          sumt2 = 0.
          do iy = iy1,iy4
            if((extwt(iy)/=0.).and.(illum(nx2,iy)==1)) then
              templw = templ(nx2,iy)*extwt(iy)
              sumt = sumt+abs(templw)
              sumt2 = sumt2+templw**2
            endif
          enddo   
          if(sumt2>0.) then 
            wtnorm = omega*sumt/sumt2
          else
            wtnorm = 0.
          endif
          do iy = iy1,iy4
            wtfac = extwt(iy)*wtnorm
            do ix = 1,nx
              templ(ix,iy) = templ(ix,iy)*wtfac
            enddo   
          enddo   
          do ix = 1,nx
            extsum = 0.
            xllsum = 0.
            do iy = iy1,iy4
              if(extwt(iy)==0.) then
                extsum = extsum+1.
                if(illum(ix,iy)==1) xllsum = xllsum+1.
              endif
            enddo
            if(xllsum<extsum) then
              if(xllsum<0.75*extsum) then
                renorm = 0.
              else
                renorm = extsum/xllsum
              endif
              do iy = iy1,iy4
                templ(ix,iy) = renorm*templ(ix,iy)
              enddo
            endif
          enddo
        else
          do iy = iy1,iy4
            templi = omega*extwt(iy)
            do ix = 1,nx
              if(illum(ix,iy)>=1) then
                templ(ix,iy) = templi
              else
                templ(ix,iy) = 0.
              endif
            enddo   
          enddo   
          do ix = 1,nx
            extsum = 0.
            xllsum = 0.
            do iy = iy1,iy4
              if(extwt(iy)==0.) then
                extsum = extsum+1.
                if(illum(ix,iy)==1) xllsum = xllsum+1.
              endif
            enddo
            if(xllsum<extsum) then
              if(xllsum<0.75*extsum) then
                renorm = 0.
              else
                renorm = extsum/xllsum
              endif
              do iy = iy1,iy4
                templ(ix,iy) = renorm*templ(ix,iy)
              enddo
            endif
          enddo
        endif

!	weight and sum spatially

  380   dlnw = pixelwd/(2.*xdr*xdfl)
        specmax = 0.
        atmax = 0.
        do ix = 1,nx
          x = ix-nx/2
          sum = 0.
          sumvar = 0.
          sumcard = 0.
          sumblk = 0.
          sumwt = 0.
          do iy = iy1,iy4
            templi = templ(ix,iy)
            templa = abs(templi)
            if(.not.justwno) then
              sum = sum+templi*arr(ix,iy)
              if(scan) then
                sumvar = sumvar+templa*std(ix,iy)
              else
                sumvar = sumvar+(templi*std(ix,iy))**2
              endif
            endif
            if((modeobs==MSTARE).and.(modecard==MNONE)) then
              sumcard = sumcard+templa*arr(ix,iy)
              sumblk = sumblk+templa
            else
              sumcard = sumcard+templa*cards(ix,iy,1)
              sumblk = sumblk+templa*cards(ix,iy,3)
            endif
            sumwt = sumwt+templa
          enddo   
          arr(ix,1) = wnoc*sexp(dlnw*x)
          if(sumwt==0.) then
            arr(ix,2) = 0.
            arr(ix,3) = 0.
            arr(ix,4) = 0.
            black(ix) = 0.
          else
            arr(ix,2) = sum
            if(scan) then
              arr(ix,3) = sumvar/sumwt
            else
              arr(ix,3) = sqrt(sumvar)
            endif
            arr(ix,4) = sumcard/sumwt
            black(ix) = sumblk/sumwt
            if(sum>specmax) specmax  = sum
            if(arr(ix,4)>atmax) atmax = arr(ix,4)
          endif
        enddo   

        if(modecard==MNONE) goto 399

        nplt = 1
        scale = 1.
        if(wno0==0.) then
          nsmoo = 2
          smax = 0.
  388     do ix = 1,nx
              plot(ix,1) = arr(ix,1)
              plot(ix,2) = arr(ix,4)
              plot(ix,3) = 0.
              plot(ix,4) = 0.
              if(arr(ix,4)>smax) smax = arr(ix,4)
          enddo   
          scale = smax
          if(ask) scale = 0.
          print '(" Plotting atmospheric transmission")'
          call plots(ms,nx,nplt,plot,scale,nomon,iwin2)
          print '(" Enter reference wavenumber ",$)'
          read *, wnoref
          if(wnoref<=0.) goto 399
          print '(" Click cursor on reference line")'
          call getcurs(x,y,ikey,ierr)
!	 w = x
          if((modecard==MSKY).or.(modecard==MNONE)) then
            w = plmax(plot,ms,x,nsmoo)
          else
            w = plmin(plot,ms,x,nsmoo)
          endif
          wnoc = wnoc+wnoref-w
          wno0 = wnoc
          print '(" wno0 = ",f10.3,"  OK? (y) ",$)', wno0
          read '(a8)', yn
          sinang = lorder/(2.0*xddgr*wno0)
          xdr = sinang/sqrt(1.-sinang**2)
          dlnw = pixelwd/(2.*xdr*xdfl)
          do ix = 1,nx
            x = ix-nx/2
            arr(ix,1) = wnoc*sexp(dlnw*x)
          enddo   
          if(index(yn,'n')>0) goto 388
          write(iupl,'("wno0 = ",f10.3)') wno0
          if(scan) ierr = 1
        elseif((wno0<0.).and.(modecard/=MNONE)) then
!       read atmo file and interpolate onto spec wno grid
          wno0 = abs(wno0)
          if(wno0<500) then
            atmofile = 'nocando'
          elseif(wno0<700) then
            atmofile = atmodir(1:lenatmo)//'Q.mod'
          elseif(wno0<1000) then
            atmofile = atmodir(1:lenatmo)//'NL.mod'
          elseif(wno0<1400.) then
            atmofile = atmodir(1:lenatmo)//'NS.mod'
          elseif(wno0<1800.) then
            atmofile = 'nocando'
          elseif(wno0<2300.) then
            atmofile = atmodir(1:lenatmo)//'M.mod'
          else
            atmofile = 'nocando'
          endif
          open(unit=iuatmo,file=atmofile,access='sequential',iostat=jerr)
          if(jerr/=0) then
            print '(" Error opening ",a32)', atmofile
            goto 399
          endif
          dwnoat = 0.001
          dwnoar = (arr(nx,1)-arr(1,1))/(nx-1)
          airpow = max(1.,min(2.,airmass))
          wnoix = arr(1,1)-dwnoar/2.
          do
            read(unit=iuatmo,fmt=*,iostat=jerr) &
              wnoa,trhrlw,trhrhw,bshrlw,bshrhw,bslrlw,bslrhw
            if(jerr/=0) then
              print '(" Error reading atmo model",2f8.3)', wnoix,wnoa
              write(iupl,'("!!! Error reading atmo model",2f8.3)') &
                wnoix,wnoa
              goto 399
            endif
            if(wnoa>wnoix-dwnoat) exit
          enddo
          do ix = 1,nx
            wnoix = arr(ix,1)-dwnoar/2.
            wnojx = wnoix+dwnoar
            sumbslw = 0.
            sumbshw = 0.
            sumn = 0.
            do
              read(unit=iuatmo,fmt=*,iostat=jerr) &
                wnoa,trhrlw,trhrhw,bshrlw,bshrhw,bslrlw,bslrhw
              if(jerr/=0) then
                print '(" Error reading atmo model",2f8.3)', wnoix,wnoa
                write(iupl,'("!!! Error reading atmo model",2f8.3)') &
                  wnoix,wnoa
                goto 399
              endif
              sumbslw = sumbslw+bslrlw
              sumbshw = sumbshw+bslrhw
              sumn = sumn+1.
              if(wnoa>wnojx-dwnoat) exit
            enddo
            bslow = sumbslw/sumn
            bshiw = sumbshw/sumn
            plot(ix,1) = arr(ix,1)
            plot(ix,2) = arr(ix,4)
            plot(ix,3) = bshiw**airpow
            plot(ix,4) = (bshiw*bslow)**(airpow/2.)
            plot(ix,5) = bslow**airpow
          enddo
          close(unit=iuatmo)
          print '(" Plotting blk-sky, models")'
          call plot4(ms,nx,nplt,plot,scale,nomon,iwin2)
          if(testrun) then
            call waitasec(.true.,.false.,plot,ms,ms,nx)
          else
            call waitasec(pause,.false.,plot,ms,ms,nx)
          endif
!	subtract smoothed spectra
          do ipl = 2,5
            do ix = 1,nx
              kx = min(16,ix-1,nx-ix)
              xk = kx+1.
              sumpl = 0.
              sumwt = 0.
              do jdx = -kx,kx
                jx = ix+jdx
                if(plot(jx,2)/=0.) then
                  wtj = 1.-abs(jdx)/xk
                  sumpl = sumpl+wtj*plot(jx,ipl)
                  sumwt = sumwt+wtj
                endif
              enddo
              if(sumwt>0.) then
                plsm(ix) = sumpl/sumwt
              else
                plsm(ix) = 0.
              endif
            enddo
            do ix = 1,nx
              if(plot(ix,2)/=0.) then
                plot(ix,ipl) = plot(ix,ipl)-plsm(ix)
              else
                plot(ix,ipl) = 0.
              endif
            enddo
          enddo
!	find shift between blk-sky and atmo model
          idxmin = -8
          chisqmin = 1.0e+06
          do idx = -7,7
            sumerr3 = 0.
            sumsq3 = 0.
            sumerr4 = 0.
            sumsq4 = 0.
            sumerr5 = 0.
            sumsq5 = 0.
            sumn = 0.
            do ix = 8,nx-7
              if(plot(ix,2)/=0.) then
                dbs = plot(ix,2)-plot(ix+idx,3)
                sumerr3 = sumerr3+dbs
                sumsq3 = sumsq3+dbs**2
                dbs = plot(ix,2)-plot(ix+idx,4)
                sumerr4 = sumerr4+dbs
                sumsq4 = sumsq4+dbs**2
                dbs = plot(ix,2)-plot(ix+idx,5)
                sumerr5 = sumerr5+dbs
                sumsq5 = sumsq5+dbs**2
                sumn = sumn+1.
              endif
            enddo
            if(sumn==0.) sumn = 1.
            chisq3 = sumsq3-sumerr3**2/sumn
            chisq4 = sumsq4-sumerr4**2/sumn
            chisq5 = sumsq5-sumerr5**2/sumn
            slope = (chisq5-chisq3)/2.
            curve = chisq3-2.*chisq4+chisq5
            if(slope==0.) then
              print '(" Bad chisq for pixel shift ",i2,5f9.3)', &
                idx,chisq3,chisq4,chisq5,slope,curve
              chisq(idx+8) = 100.
              cycle
            endif
            dmin = -slope/max(abs(slope),curve)
            chisqi = chisq4+slope*dmin+curve*dmin**2/2.
            chisq(idx+8) = chisqi
            if(chisqi<chisqmin) then
              idxmin = idx
              chisqmin = chisqi
            endif
          enddo
          if(abs(idxmin)>6) then
            print '(" Best pixel shift > 6")'
            if(doplot) then
              print '(" Switch to kbd mode? (n) "$)'
              read '(a8)', yn
              if(index(yn,'y')>0) then
                ierr = 9
                return
              endif
            endif
          else
            chisq1 = chisq(idxmin+7)
            chisq2 = chisq(idxmin+8)
            chisq3 = chisq(idxmin+9)
            didx = idxmin+((chisq1-chisq3)/2.)/(chisq1+chisq3-2.*chisq2)
            print '(" Best pixel shift =",f6.2)', didx
            if(abs(didx-idxmin)>1.) &
              print '(8f10.3/7f10.3,i6)', (chisq(idx),idx=1,15),idxmin
          endif
          if(abs(idxmin)>2) then
  385       print '(" Shift by p(ixel), w(no), n(one), or s(kip) "$)'
            read '(a8)', yn
            if(index(yn,'n')>0) then
              didx = 0.
              goto 390
            elseif(index(yn,'s')>0) then
              ierr = 8
              return
            else
              print '(" Enter shift ", & 
                "(positive to shift observed spectrum blue) ",$)'
              read(*,*,iostat=ierr) shift
              if(ierr>0) goto 385
              if(shift==0.) then
                didx = 0.
                goto 390
              elseif(index(yn,'p')>0) then
                wno0 = wno0+shift*dwnoar
              elseif(index(yn,'w')>0) then
                wno0 = wno0+shift
              else
                goto 385
              endif
              wnoc = wno0
              wno0 = -wno0
              goto 380
            endif
          endif
  390     wno0 = abs(wno0)+didx*dwnoar
          wnoc = wno0
          print '(" Pixel shift, new wno0:",2f10.4)', didx,wno0
          write(iupl,'(" Pixel shift, new wno0:",f10.4)') didx,wno0
          sinang = lorder/(2.0*xddgr*wno0)
          if((abs(sinang)>0.99).or.(0.*sinang/=0.)) then
            print '(" order,xddgr,wno0,sinang =",i6,4es10.2)', &
              lorder,xddgr,wno0,sinang
            print '(" Enter sin(grating) ",$)'
            read *, sinang
          endif
          xdr = sinang/sqrt(1.-sinang**2)
          dlnw = pixelwd/(2.*xdr*xdfl)
          do ix = 1,nx
            x = ix-nx2
            arr(ix,1) = wno0*sexp(dlnw*x)
          enddo
          if(modeinst==MMED) then
            xdang = echelle
          else
            xdang = lores
          endif
          xdorder = nint(2.0*xddgr*wno0*sin(xdang/DEGRAD))
          xdang0 = DEGRAD*asin(xdorder/(2.0*xddgr*wno0))
          if(modeinst==MMED) then
            xdmr0 = echelle+xdmr0-xdang0
            echelle = xdang0
            print '(" mr0 = ",f10.3)', xdmr0
            write(iupl,'(" mr0 = ",f10.3)') xdmr0
          else
            xdlr0 = lores+xdlr0-xdang0
            lores = xdang0
            print '(" lr0 = ",f10.3)', xdlr0
            write(iupl,'(" lr0 = ",f10.3)') xdlr0
          endif
          write(iuwno,'(a24,2f8.3)') rawfile,wno0
        else
          do ix = 1,nx
            plot(ix,1) = arr(ix,1)
            plot(ix,2) = arr(ix,4)
            plot(ix,3) = 0.
            plot(ix,4) = 0.
            plot(ix,5) = 0.
          enddo
          print '(" Plotting reference spectrum")'
          call plot4(ms,nx,nplt,plot,scale,nomon,iwin2)
          call waitasec(pause,.false.,plot,ms,ms,nx)
        endif
        hrr = 9.8

  399   if(justwno.or.(modeext==MNONE)) goto 369
        sumsn = 0.
        sumn = 0.
        if(specmax==0.) specmax = 1.
        if((modecard/=MBLKSKY).and.(modecard/=MBSBS)) &
          specmax = specmax/atmax
        dopp = (1.+radvel/CLIGHT)
        do ix = 1,nx
          plot(ix,1) = dopp*arr(ix,1)
          plot(ix,2) = arr(ix,2)
          plot(ix,3) = arr(ix,3)
          plot(ix,4) = arr(ix,4)
          plot(ix,5) = black(ix)
          if((arr(ix,2)*arr(ix,3))/=0.) then
            sumsn = sumsn+arr(ix,2)/arr(ix,3)
            sumn = sumn+1.
          endif
        enddo   
        sumsn = sumsn/sumn
        print '(" Mean S/N in extracted spectrum: ",f10.2)', sumsn
        write(iupl,'(" Mean S/N in extracted spectrum: ",f10.2)') sumsn

        scale = -1.
        if(ask) scale = 0.
        nplt = 1
        print '(" Plotting extracted spectrum")'
        if(domon) then
          call plots(ms,nx,nplt,plot,scale,outfile,iwin2)
        else
          call plots(ms,nx,nplt,plot,scale,nomon,iwin2)
        endif
        call waitasec(pause,.false.,plot,ms,ms,nx)
  369   continue

        wnomin = arr(1,1)
        wnomax = arr(nx,1)
        dvpix = -CLIGHT*dlnw
        rslit = pltscl/(dlnw*slitwid)
        rdiff = 13.5*xdr*wno0
        resolv = 1.0/sqrt(1.0/rslit**2+1.0/rdiff**2)
        inquire(unit=iuredh,opened=hdopen)
        if(hdopen) then
          write(iuredh,'("wno0    = ",f10.4)') wnoc
          if(dofits) then
            comment = 'wavenumber at pixel 128'
            call fithreal('WNO0    ',wno0,comment,iufith)
            comment = 'minimum wavenumber (cm-1)'
            call fithreal('WNO_MIN ',wnomin,comment,iufith)
            comment = 'maximum wavenumber (cm-1)'
            call fithreal('WNO_MAX ',wnomax,comment,iufith)
            comment = 'velocity interval between pixels (km/s)'
            call fithreal('DVPIX   ',dvpix,comment,iufith)
            comment = 'approximate resolving power'
            call fithreal('RESOLV  ',resolv,comment,iufith)
          endif
        endif
        inquire(unit=iusumh,opened=hdopen)
        if(dosum.and.(sumtime==0.).and.hdopen)  then
          write(iusumh,'("wno0    = ",f10.4)') wnoc
          if(dofits) then
            comment = 'wavenumber at pixel 128'
            call fithreal('WNO0    ',wno0,comment,iufish)
            comment = 'minimum wavenumber (cm-1)'
            call fithreal('WNO_MIN ',wnomin,comment,iufish)
            comment = 'maximum wavenumber (cm-1)'
            call fithreal('WNO_MAX ',wnomax,comment,iufish)
            comment = 'velocity interval between pixels (km/s)'
            call fithreal('DVPIX   ',dvpix,comment,iufish)
            comment = 'approximate resolving power'
            call fithreal('RESOLV  ',resolv,comment,iufish)
          endif
        endif

      endif

!        if(justwno) then
!          close(unit=iufith)
!          close(unit=iuredh)
!        endif

        return

      end


      subroutine zerosum(ierr)
!	initialize specsum and arrsum to zero

        use dims

        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime

        do it = 1,mt
          do is = 1,ms
            specsum(is,it) = 0.
          enddo   
        enddo   

        do iy = 1,my
          do ix = 1,mx
            arrsum(ix,iy) = 0.
          enddo   
        enddo   

        do ip = 1,mp
          do iy = 1,my
            do ix = 1,mx
              scansum(ix,iy,ip) = 0.
            enddo   
          enddo   
        enddo   

        sumtime = 0.

        ierr = 0
        return
      end


      subroutine sumarr(arr,ierr)
!	average long-slit array with arrsum

        use dims

        real arr(mx,my)
        real nodpa,lores,kmirror,krot

        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl

        oldtime = sumtime
        sumtime = sumtime+addtime
        t1 = oldtime/sumtime
        t2 = addtime/sumtime
        do iy = 1,my
          do ix = 1,mx
            if((arrsum(ix,iy)/=0.).or.(t1==0.)) then
              if(arr(ix,iy)/=0.) then
                arrsum(ix,iy) = t1*arrsum(ix,iy)+t2*arr(ix,iy)
              else
                arrsum(ix,iy) = 0.
              endif
            endif
          enddo   
        enddo   

        ierr = 0
        return
      end


      subroutine sumspec(spec,ierr)
!	average xd spectrum with specsum

        use dims

        real spec(ms,mt)
        real nodpa,lores,kmirror,krot

        common /nn/ nx,ny,nc,norder,ns,nt
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl

        oldtime = sumtime
        sumtime = sumtime+addtime
        t1 = oldtime/sumtime
        t2 = addtime/sumtime
        do it = 1,mt
          do is = 1,ms
            if((specsum(is,it)/=0.).or.(t1==0.)) then
              if(spec(is,it)/=0.) then
                if(it==3) then
                  specsum(is,it) = &
                    sqrt((t1*specsum(is,it))**2+(t2*spec(is,it))**2)
                else
                  specsum(is,it) = t1*specsum(is,it)+t2*spec(is,it)
                endif
              else
                specsum(is,it) = 0.
              endif
            endif
          enddo   
        enddo   

        ierr = 0
        return
      end


      subroutine sumscan(scan,scansum,std,nz,im,sumtime,ierr)
!	clean up a scan, optionally shift, and add to scansum

        use ius
        use dims
        use modes

        real scan(mx,my,1),scansum(mx,my,1),std(mx,my)
        real corrx(mp),corrz(mp),sky(mx,mp)
        integer jllx(mx)
        real nodpa,lores,kmirror,krot
        logical baddata,doplot
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical contwt,badx,savask
        character(16) yn
        character(60) comment

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /iwin/ iwin1,iwin2

        save ix1,ix2,iz1,iz2,ip1,ip2,it1,it2,iorder1,iorder2

        doplot = (iwin1>0)
        nzscan = nz-nsky-1
        badx = .false.

!	set scan(ix,iy,iz) = 0 where illum != 1
!	and jllx(ix) = 0 outside of shiftint

        if(crossdisp) then
          nsh = min0(2,nt/3)
        else
          nsh = 4
        endif
        nsh2 = nsh/2
!	it appears that the following is overwritten below
        if(im==1) then
          if(crossdisp) then
            if((intshift(1)*intshift(2))>0) then
              it1 = max(1,abs(intshift(1)))+nsh
              it2 = min(nt,abs(intshift(2)))-nsh
            else
              it1 = 1+nsh
              it2 = nt-nsh
            endif
            if((intshift(3)*intshift(4))>0) then
              iorder1 = max(1,intshift(3))
              iorder2 = min(norder,intshift(4))
            else
              iorder1 = 1
              iorder2 = norder
            endif
            ix1 = it1+(iorder1-1)*nt
            ix2 = it2+(iorder2-1)*nt
          else
            if((intshift(1)*intshift(2))>0) then
              it1 = max(1,intshift(1))+nsh
              it2 = min(ny,intshift(2))-nsh
            else
              it1 = 1+nsh
              it2 = ny-nsh
            endif
          endif
          if((intshift(5)*intshift(6))>0) then
            ip1 = max(1,abs(intshift(5)))
            ip2 = min(ny,abs(intshift(6)))
            contwt = (intshift(5)>0)
          else
            ip1 = 1
            ip2 = ny
            contwt = .true.
          endif
          if((intsky(1)*intsky(2))>0) then
            iz1 = intsky(2)+nsh
            if((intsky(3)*intsky(4))>0) then
              iz2 = intsky(3)-nsh
            else
              iz2 = nzscan-nsh
            endif
          else
            iz1 = 3+nsh
            iz2 = nzscan-nsh
          endif
          if(doshift.and.((iz1+3)>iz2)) then
            if(nzscan>(2*nsh+3)) then
              iz1 = 3+nsh
              iz2 = nzscan-nsh
            else
              print '("*** Scan too short to shift")'
              doshift = .false.
            endif
          endif
        endif

        njll = 0
        do ix = 1,nx
          do iy = 1,ny
            if(illum(ix,iy)<=0) then
              do iz = 1,nz
                scan(ix,iy,iz) = 0.
              enddo   
              std(ix,iy) = 0.
            endif
          enddo   
          if(crossdisp) then
            jllx(ix) = 0
            it = mod(ix-1,nt)+1
            iorder = (ix-1)/nt+1
            if((it<it1).or.(it>it2) &
              .or.(iorder<iorder1).or.(iorder>iorder2) &
              .or.((ix-nsh)<1).or.((ix+nsh)>nx)) cycle
            do ish = -nsh2,nsh2
              if(illx(ix+ish)<=0) cycle
            enddo   
            jllx(ix) = 1
            njll = njll+1
          else
            jllx(ix) = 0
            if(((ix<it1).or.(ix>it2)) &
              .or.((ix-nsh)<1).or.((ix+nsh)>ny)) cycle
            do ish = -nsh,nsh
              if(illx(ix+ish)<=0) cycle
            enddo   
            jllx(ix) = 1
            njll = njll+1
          endif
        enddo   
        if(doshift.and.(njll==0)) then
          print '("*** Error in scan correlation.", &
            " No illuminated pixels in shiftint.")'
          print '(5i4)', it1,it2,nsh,iorder1,iorder2
          print '(36i2)', (illx(ix),ix=1,nt)
        endif

!	make postage stamp images and find peak

      if(doplot.and.(intshift(1)>=0)) then
        savask = ask
        if(doshift.and.(im==1)) then
          print '(" Plotting cleaned scan map")'
          print '(" Find best peak and width in order",2i4)', &
            iorder1,iorder2
          ask = .true.
        elseif(doshift.and.ask) then
          print '(" Find shifts")'
        endif
        nzoff = nzscan
        call scanmap2(scan,scansum,std,nzoff)
        ask = savask
        print '(" Sumscan is offset by",i4," pixels")', nzoff
        if(ask.and.(.not.doshift)) then
          print '(" Keep scan? (y) ",$)'
          read '(a8)', yn
          if(index(yn,'n')>0) return
        elseif(doshift.and.doplot.and.(im==1)) then
          print '(" Enter peak x,z, FWHM x,z")'
          read *, pkx,pkz,fwx,fwz
          if(pkx*pkz*fwx*fwz/=0.) then
            jx1 = nint(pkx-fwx)
            jx2 = nint(pkx+fwx)
            jorder = (jx1+jx2-1)/(2*nt)+1
            it1 = max(it1,jx1-(jorder-1)*nt)
            it2 = min(it2,jx2-(jorder-1)*nt)
            if(jorder>norder) then
              ix1 = it1+(iorder1-1)*nt
              ix2 = it2+(iorder2-1)*nt
            else
              ix1 = max(ix1,jx1)
              ix2 = min(ix2,jx2)
            endif
            iz1 = max(iz1,nint(pkz-fwz))
            iz2 = min(iz2,nint(pkz+fwz))
          endif
!          print '(8i4)', ix1,ix2,iz1,iz2,jx1,jx2,nz,nsky
        endif
      endif

!	correlate scan with scansum and shift

      if(crossdisp) then

!	x correlation

        ishx = 0
        ishz = 0
        if(doshift.and.(im>1)) then
          nsh = min0(4,nt/3)
          if(nsh<4) print '(" Allowing max xshift of",i2)', nsh-1
          xnsh = nsh-0.5
          iy1 = ip1
          iy2 = ip2

          imx = 0
          corrmx = 0.
          do ish = -nsh,nsh
            corr = 0.
            do iz = iz1,iz2
              do ix = ix1,ix2
                if(jllx(ix)/=1) cycle
                ixsh = ix-ish
                if(contwt) then
                  scancont = 0.
                  sumcont = 0.
                  do iy = iy1,iy2
                    if(illum(ix,iy)==1) then
                      if(illum(ixsh,iy)==1) &
                        scancont = scancont+scan(ixsh,iy,iz)/std(ixsh,iy)
                      sumcont = sumcont+scansum(ix,iy,iz)/std(ix,iy)
                    endif
                  enddo   
!	          corr = corr+scancont*sumcont
                  corr = corr-(scancont-sumcont)**2
                else
                  do iy = iy1,iy2
                    if(illum(ix,iy)==1) then
                      if(illum(ixsh,iy)==1) then
                        corr = corr &
                  -((scan(ixsh,iy,iz)-scansum(ix,iy,iz))/std(ix,iy))**2
                      else
                        corr = corr-(scansum(ix,iy,iz)/std(ix,iy))**2
                      endif
                    endif
                  enddo   
                endif
              enddo   
            enddo   
            corrx(ish+nsh+1) = corr
            if((corr>corrmx).or.(corrmx==0.)) then
              corrmx = corr
              imx = ish
            endif
          enddo   
          if(corrmx==0.) then
            print '("*** Error in scan correlation.")'
            print '(6i4)', ix1,ix2,iy1,iy2,iz1,iz2
            print '(8es10.2)', (corrx(ish),ish=1,2*nsh)
          endif

!	find peak of correlation

          imx = max(-nsh+1,min(nsh-1,imx))
          p1 = corrx(imx+nsh)
          p2 = corrx(imx+nsh+1)
          p3 = corrx(imx+nsh+2)
          pa = p2
          pb = (p3-p1)/2.
          pc = (p1+p3)/2.-p2
          xmax = imx-pb/(2.*pc)
          print '(" Calculated x shift = ",f6.2)', xmax
          if((0.*xmax)/=0.) print *, imx,corrmx,p1,p2,p3,pa,pb,pc
          if((corrmx==0.).or.(.not.(abs(xmax)<=xnsh))) then
            if(doplot) then
              print '(9f8.3)', (corrx(i)/abs(corrmx),i=1,9)
              print '(" Enter x shift or 99 to skip scan ",$)'
              read *, xmax
            else
              print '("*** skipping scan")'
              xmax = 99.
            endif
            badx = .true.
          elseif(ask) then
            print '(" OK?  (y) ",$)'
            read '(a16)', yn
            if(index(yn,'n')>0) then
              print '(" Enter x shift or 99 to skip scan ",$)'
              read *, xmax
            else
              read(unit=yn,fmt=*,err=109,end=109) xmax
            endif
          endif
  109     if(abs(xmax)>=99.) then
            im = im-1
            if(dosum) then
              write(iusumh,'("skip    = ",i4)') im
              comment = ' '
              if(dofits) call fithint('SKIP    ',im,comment,iufish)
            else
              write(iuredh,'("skip    = ",i4)') im
              comment = ' '
              if(dofits) call fithint('SKIP    ',im,comment,iufith)
            endif
            return
          endif
          imx = nint(xmax)

!	shift scan in x

          ishx = imx-imx/im
          if(ishx/=0) then
            if(ishx<=0) then
              nx1 = 1
              nx2 = nx+ishx
              idx = 1
            else
              nx1 = nx
              nx2 = 1+ishx
              idx = -1
            endif
            do iz = 1,nz
              do iy = 1,ny
                do ix = nx1,nx2,idx
                  do jx = ix,ix-ishx,idx
                    if(illx(jx)<=0) then
                      scan(ix,iy,iz) = 0.
                      goto 116
                    endif
                  enddo   
                  scan(ix,iy,iz) = scan(ix-ishx,iy,iz)
  116           enddo   
                if(ishx<0) then
                  do ix = nx2+1,nx
                    scan(ix,iy,iz) = 0.
                  enddo   
                else
                  do ix = 1,nx2-1
                    scan(ix,iy,iz) = 0.
                  enddo   
                endif
              enddo   
            enddo   
          endif
          ishx = imx

!	z correlation

          nsh = min0(4,nzscan/3)
          if(nsh<4) print '(" Allowing max zshift of",i2)', nsh-1
          xnsh = nsh-0.5
          imx = 0
          corrmx = 0.
          do ish = -nsh,nsh
            corr = 0.
            do iz = iz1,iz2
              izsh = iz-ish
              do ix = ix1,ix2
                if(jllx(ix)/=1) cycle
                if(contwt) then
                  scancont = 0.
                  sumcont = 0.
                  do iy = iy1,iy2
                    if(illum(ix,iy)==1) then
                      stdi = std(ix,iy)
                      scancont = scancont+scan(ix,iy,izsh)/stdi
                      sumcont = sumcont+scansum(ix,iy,iz)/stdi
                    endif
                  enddo   
!	          corr = corr+scancont*sumcont
                  corr = corr-(scancont-sumcont)**2
                else
                  do iy = iy1,iy2
                    if(illum(ix,iy)==1) &
                      corr = corr &
                  -((scan(ix,iy,izsh)-scansum(ix,iy,iz))/std(ix,iy))**2
                  enddo   
                endif
              enddo   
            enddo   
            corrz(ish+nsh+1) = corr
            if((corr>corrmx).or.(corrmx==0.)) then
              corrmx = corr
              imx = ish
            endif
          enddo   

          imx = max(-nsh+1,min(nsh-1,imx))
          p1 = corrz(imx+nsh)
          p2 = corrz(imx+nsh+1)
          p3 = corrz(imx+nsh+2)
          pa = p2
          pb = (p3-p1)/2.
          pc = (p1+p3)/2.-p2
          xmax = imx-pb/(2.*pc)
          if((0.*xmax)/=0.) print *, imx,corrmx,p1,p2,p3,pa,pb,pc
          print '(" Calculated z shift = ",f6.2)', xmax
          if(.not.(abs(xmax)<=xnsh)) then
            if(doplot) then
              print '(9f8.3)', (corrz(i)/abs(corrmx),i=1,9)
              print '(" Enter z shift or 99 to skip scan ",$)'
              read *, xmax
            else
              xmax = 99.
            endif
          elseif(ask.or.badx) then
            print '(" OK?  (y) ",$)'
            read '(a16)', yn
            if(index(yn,'n')>0) then
              print '(" Enter z shift or 99 to skip scan ",$)'
              read *, xmax
            else
              read(unit=yn,fmt=*,err=209,end=209) xmax
            endif
          endif
  209     if(abs(xmax)>=99.) then
            im = im-1
            if(dosum) then
              write(iusumh,'("skip    = ",i4)') im
              comment = ' '
              if(dofits) call fithint('SKIP    ',im,comment,iufish)
            else
              write(iuredh,'("skip    = ",i4)') im
              comment = ' '
              if(dofits) call fithint('SKIP    ',im,comment,iufith)
            endif
            return
          endif
          imx = nint(xmax)

!	shift scan in z

          ishz = imx-(imx/im)
          if(ishz/=0) then
            if(ishz<0) then
              nz1 = 1
              nz2 = nzscan+ishz
              idz = 1
            else
              nz1 = nzscan
              nz2 = 1+ishz
              idz = -1
            endif
            do iy = 1,ny
              do ix = 1,nx
                do iz = nz1,nz2,idz
                  if(illum(ix,iy)==1) then
                    scan(ix,iy,iz) = scan(ix,iy,iz-ishz)
                  else
                    scan(ix,iy,iz) = 0.
                  endif
                enddo   
!                if(ishz<0) then
!                  do iz = nz2+1,nzscan
!                    scan(ix,iy,iz) = 0.
!                  enddo   
!                else
!                  do iz = 1,nz2-1
!                    scan(ix,iy,iz) = 0.
!                  enddo   
!                endif
              enddo   
            enddo   
          endif
          ishz = imx
          print '(" x,z shifts for scan",i3,":",2i3)', im,ishx,ishz
          write(iupl,'(" x,z shifts for scan",i3,":",2i3)') im,ishx,ishz
          if(dosum) then
            write(iusumh,'("yshift  = ",i4)') ishx
            write(iusumh,'("zshift  = ",i4)') ishz
            if(dofits) then
              write(unit=comment,fmt='("for scan",i3)') im
              call fithint('YSHIFT  ',ishx,comment,iufish)
              call fithint('ZSHIFT  ',ishz,comment,iufish)
            endif
          else
            write(iuredh,'("yshift  = ",i4)') ishx
            write(iuredh,'("zshift  = ",i4)') ishz
            if(dofits) then
              write(unit=comment,fmt='("for scan",i3)') im
              call fithint('YSHIFT  ',ishx,comment,iufith)
              call fithint('ZSHIFT  ',ishz,comment,iufith)
            endif
          endif

        endif

!	add scan to scansum, applying remaining shift to scansum

        if(im>1) then
          x2 = 1./im
          ishx = ishx/im
          ishz = ishz/im
        else
          x2 = 1.
          ishx = 0
          ishz = 0
        endif
        x1 = 1.-x2
        if(ishx>=0) then
          nx1 = 1
          nx2 = nx-ishx
          nx3 = nx2+1
          nx4 = nx
          idx = 1
        else
          nx1 = nx
          nx2 = 1-ishx
          nx3 = nx2-1
          nx4 = 1
          idx = -1
        endif
        if(ishz>=0) then
          nz1 = 1
          nz2 = nzscan-ishz
          nz3 = nz2+1
          nz4 = nzscan
          idz = 1
        else
          nz1 = nzscan
          nz2 = 1-ishz
          nz3 = nz2-1
          nz4 = 1
          idz = -1
        endif
        do iz = nz1,nz2,idz
          izsh = iz+ishz
          do iy = 1,ny
            do ix = nx1,nx2,idx
              ixsh = ix+ishx
              if((x1==0.).or.(scansum(ixsh,iy,izsh)==0.)) then
                scansum(ix,iy,iz) = scan(ix,iy,iz)
              elseif(scan(ix,iy,iz)/=0.) then
                scansum(ix,iy,iz) = &
                  x1*scansum(ixsh,iy,izsh)+x2*scan(ix,iy,iz)
              endif
            enddo   
            do ix = nx3,nx4,idx
              scansum(ix,iy,iz) = scan(ix,iy,iz)
            enddo   
          enddo   
        enddo   
        do iz = nz3,nz4,idz
          do iy = 1,ny
            do ix = 1,nx
              scansum(ix,iy,iz) = scan(ix,iy,iz)
            enddo   
          enddo   
        enddo   
        do iz = nzscan+1,nz
          do iy = 1,ny
            do ix = nx1,nx4,idx
              ixsh = ix+ishx
              if((ixsh>=1).and.(ixsh<=nx) &
                .and.(scansum(ixsh,iy,iz)/=0.)) then
                if(scan(ix,iy,iz)/=0.) &
                  scansum(ix,iy,iz) = &
                    x1*scansum(ixsh,iy,iz)+x2*scan(ix,iy,iz)
              else
                scansum(ix,iy,iz) = scan(ix,iy,iz)
              endif
            enddo   
          enddo   
        enddo   
        if(nz<=mp) then
          nz = nz+1
          iz = nz
          do iy = 1,ny
            do ix = nx1,nx4,idx
              ixsh = ix+ishx
              if((ixsh>=1).and.(ixsh<=nx) &
                .and.(scansum(ixsh,iy,iz)/=0.)) then
                if(std(ix,iy)/=0.) &
                  scansum(ix,iy,iz) = sqrt( &
                    (x1*scansum(ixsh,iy,iz))**2+(x2*std(ix,iy))**2)
              else
                scansum(ix,iy,iz) = std(ix,iy)
              endif
              scan(ix,iy,iz) = std(ix,iy)
            enddo   
          enddo   
        endif

      else
!	long-slit data

!	y correlation

        ishy = 0
        ishz = 0
        if(doshift.and.(im>1)) then
          ix1 = ip1
          ix2 = ip2
          iy1 = it1
          iy2 = it2
          nsh = 4
          xnsh = nsh-0.5
          imx = 0
          corrmx = 0.
          do ish = -nsh,nsh
            corr = 0.
            do iz = iz1,iz2
              do iy = iy1,iy2
                iysh = iy-ish
                if((jllx(iy)/=1).or.(jllx(iysh)==0)) then
                  continue
                elseif(contwt) then
                  scancont = 0.
                  sumcont = 0.
                  do ix = ix1,ix2
                    if(illum(ix,iysh)==1) &
                      scancont = scancont &
                        +scan(ix,iysh,iz)/std(ix,iysh)
                    if(illum(ix,iy)==1) &
                      sumcont = sumcont+scansum(ix,iy,iz)/std(ix,iy)
                  enddo   
!	          corr = corr+scancont*sumcont
                  corr = corr-(scancont-sumcont)**2
                else
                  do ix = ix1,ix2
                    if((jllx(iy)==1).and.(illum(ix,iy)==1)) &
                      corr = corr &
!                    +scan(ix,iysh,iz)*scansum(ix,iy,iz)/std(ix,iy)**2 &
                      -(scan(ix,iysh,iz)-scansum(ix,iy,iz))**2 &
                      /(std(ix,iy)*std(ix,iysh))
                  enddo   
                endif
              enddo   
            enddo   
            corrx(ish+nsh+1) = corr
            if((corr>corrmx).or.(corrmx==0.)) then
              corrmx = corr
              imx = ish
            endif
          enddo   
          if(corrmx==0.) print '("*** Error in scan correlation.")'

!	find peak of correlation

          imx = max(-nsh+1,min(nsh-1,imx))
          p1 = corrx(imx+nsh)
          p2 = corrx(imx+nsh+1)
          p3 = corrx(imx+nsh+2)
          pa = p2
          pb = (p3-p1)/2.
          pc = (p1+p3)/2.-p2
          xmax = imx-pb/(2.*pc)
          if(ibad(xmax)/=0) print *, imx,corrmx,p1,p2,p3,pa,pb,pc
          print '(" Calculated x shift = ",f6.2)', xmax
          if(.not.(abs(xmax)<=xnsh)) then
            print '(9f8.3)', (corrx(i)/abs(corrmx),i=1,9)
            print '(" Enter x shift or 99 to skip scan ",$)'
            read *, xmax
            badx = .true.
          elseif(ask) then
            print '(" OK?  (y) ",$)'
            read '(a16)', yn
            if(index(yn,'n')>0) then
              print '(" Enter x shift or 99 to skip scan ",$)'
              read *, xmax
            else
              read(unit=yn,fmt=*,err=309,end=309) xmax
            endif
          endif
  309     if(abs(xmax)>=99.) then
            im = im-1
            if(dosum) then
              write(iusumh,'("skip    = ",i4)') im
              comment = ' '
              if(dofits) call fithint('SKIP    ',im,comment,iufish)
            else
              write(iuredh,'("skip    = ",i4)') im
              comment = ' '
              if(dofits) call fithint('SKIP    ',im,comment,iufith)
            endif
            return
          endif
          imx = nint(xmax)

!	get sky from average along slit before shifting in y
          do iz = 1,nz
            do ix = 1,nx
              skysum = 0.
              illsum = 0
              do iy = iy1,iy2
                if(jllx(iy)/=0) then
                  skysum = skysum+scan(ix,iy,iz)
                  illsum = illsum+1
                endif
              enddo   
              if(illsum==0) illsum = 1
              sky(ix,iz) = skysum/illsum
            enddo   
          enddo   

!	shift scan in y, saving imx/im shift for scansum

          ishy = imx-imx/im
          if(ishy/=0) then
            if(ishy<=0) then
              ny1 = 1
              ny2 = ny+ishy
              idy = 1
            else
              ny1 = ny
              ny2 = 1+ishy
              idy = -1
            endif
            do iz = 1,nz
              do ix = 1,nx
                do iy = ny1,ny2,idy
                  do jy = iy,iy-ishy,idy
                    if(illx(jy)<=0) then
                      scan(ix,iy,iz) = 0.
                      cycle
                    endif
                  enddo   
                  scan(ix,iy,iz) = scan(ix,iy-ishy,iz)
                enddo   
                if(ishy<0) then
                  do iy = ny2+1,ny
                    scan(ix,iy,iz) = 0.
                  enddo   
                elseif(ishy>0) then
                  do iy = 1,ny2-1
                    scan(ix,iy,iz) = 0.
                  enddo   
                endif
              enddo   
            enddo   
          endif
          iy1 = max(1,it1+ishy)
          iy2 = min(ny,it2+ishy)
          ishy = imx

!	z correlation

          nsh = min0(4,nzscan/3)
          if(nsh<4) print '(" Allowing max z shift of",i2)', nsh-1
          xnsh = nsh-0.5
          imx = 0
          corrmx = 0.
          do ish = -nsh,nsh
            corr = 0.
            do iz = iz1,iz2
              izsh = iz-ish
              do iy = iy1,iy2
                if(jllx(iy)==0) then
                  continue
                elseif(contwt) then
                  scancont = 0.
                  sumcont = 0.
                  do ix = ix1,ix2
!	note: scan is shifted in y but illum and std aren't
                    if(illum(ix,iy)==1) then
                      stdi = std(ix,iy)
                      scancont = scancont &
                        +(scan(ix,iy,izsh)-sky(ix,izsh))/stdi
                      sumcont = sumcont+scansum(ix,iy,iz)/stdi
                    endif
                  enddo   
!	          corr = corr+scancont*sumcont
                  corr = corr-(scancont-sumcont)**2
                else
                  do ix = ix1,ix2
                    if(illum(ix,iy)==1) &
!                      corr = corr+(scan(ix,iy,izsh)-sky(ix,izsh)) &
!                        *scansum(ix,iy,iz)/std(ix,iy)**2
                      corr = corr-(((scan(ix,iy,izsh)-sky(ix,izsh)) &
                        -scansum(ix,iy,iz))/std(ix,iy))**2
                  enddo   
                endif
              enddo   
            enddo   
            corrz(ish+nsh+1) = corr
            if((corr>corrmx).or.(corrmx==0.)) then
              corrmx = corr
              imx = ish
            endif
          enddo   

          imx = max(-nsh+1,min(nsh-1,imx))
          p1 = corrz(imx+nsh)
          p2 = corrz(imx+nsh+1)
          p3 = corrz(imx+nsh+2)
          pa = p2
          pb = (p3-p1)/2.
          pc = (p1+p3)/2.-p2
          xmax = imx-pb/(2.*pc)
          if((0.*xmax)/=0.) print *, imx,corrmx,p1,p2,p3,pa,pb,pc
          print '(" Calculated z shift = ",f6.2)', xmax
          if(.not.(abs(xmax)<=xnsh)) then
            print '(9f8.3)', (corrz(i)/abs(corrmx),i=1,9)
            print '(" Enter z shift or 99 to skip scan ",$)'
            read *, xmax
          elseif(ask.or.badx) then
            print '(" OK?  (y) ",$)'
            read '(a16)', yn
            if(index(yn,'n')>0) then
              print '(" Enter z shift or 99 to skip scan ",$)'
              read *, xmax
            else
              read(unit=yn,fmt=*,err=409,end=409) xmax
            endif
          endif
  409     if(abs(xmax)>=99.) then
            im = im-1
            if(dosum) then
              write(iusumh,'("skip    = ",i4)') im
              comment = ' '
              if(dofits) call fithint('SKIP    ',im,comment,iufish)
            else
              write(iuredh,'("skip    = ",i4)') im
              comment = ' '
              if(dofits) call fithint('SKIP    ',im,comment,iufith)
            endif
            return
          endif
          imx = nint(xmax)

!	shift scan in z

          ishz = imx-imx/im
          if(ishz/=0) then
            if(ishz<0) then
              iz1 = 1
              iz2 = nzscan+ishz
              idz = 1
            else
              iz1 = nzscan
              iz2 = 1+ishz
              idz = -1
            endif
            do iy = 1,ny
              do ix = 1,nx
                do iz = iz1,iz2,idz
                  if(illum(ix,iy)==1) then
                    scan(ix,iy,iz) = scan(ix,iy,iz-ishz)
                  else
                    scan(ix,iy,iz) = 0.
                  endif
                enddo   
                if(ishz<0) then
                  do iz = iz2+1,nzscan
                    scan(ix,iy,iz) = 0.
                  enddo   
                else
                  do iz = 1,iz2-1
                    scan(ix,iy,iz) = 0.
                  enddo   
                endif
              enddo   
            enddo   
          endif
          ishz = imx
          print '(" y,z shifts for scan",i3,":",2i3)', im,ishy,ishz
          write(iupl,'(" y,z shifts for scan",i3,":",2i3)') im,ishy,ishz
          if(dosum) then
            write(iusumh,'("yshift  = ",i4)') ishx
            write(iusumh,'("zshift  = ",i4)') ishz
            if(dofits) then
              write(unit=comment,fmt='("for scan",i3)') im
              call fithint('YSHIFT  ',ishx,comment,iufish)
              call fithint('ZSHIFT  ',ishz,comment,iufish)
            endif
          else
            write(iuredh,'("yshift  = ",i4)') ishx
            write(iuredh,'("zshift  = ",i4)') ishz
            if(dofits) then
              write(unit=comment,fmt='("for scan",i3)') im
              call fithint('YSHIFT  ',ishx,comment,iufith)
              call fithint('ZSHIFT  ',ishz,comment,iufith)
            endif
          endif

        endif

!	add scan to scansum, applying remaining shift to scansum

        if(im>1) then
          xim = im
          x20 = 1./xim
          x21 = 1./(xim-1.)
          ishy = ishy/im
          ishz = ishz/im
        else
          x20 = 1.
          x21 = 1.
          ishy = 0
          ishz = 0
        endif
        x10 = 1.-x20
        x11 = 1.-x21
        if(ishy>=0) then
          iy1 = 1
          iy2 = ny-ishy
          iy3 = iy2+1
          iy4 = ny
          idy = 1
        else
          iy1 = ny
          iy2 = 1-ishy
          iy3 = iy2-1
          iy4 = 1
          idy = -1
        endif
        if(ishz>=0) then
          iz1 = 1
          iz2 = nzscan-ishz
          iz3 = iz2+1
          iz4 = nzscan
          idz = 1
        else
          iz1 = nzscan
          iz2 = 1-ishz
          iz3 = iz2-1
          iz4 = 1
          idz = -1
        endif
        do iz = iz1,iz2,idz
!	to avoid the memory of the card
!	  if(iz<=2) then
!	    x1 = x11
!	    x2 = x21
!	  else
            x1 = x10
            x2 = x20
!	  endif
          izsh = iz+ishz
          do iy = iy1,iy2,idy
            iysh = iy+ishy
            do ix = 1,nx
              if((x1==0.).or.(scansum(ix,iysh,izsh)==0.)) then
                scansum(ix,iy,iz) = scan(ix,iy,iz)
              elseif(scan(ix,iy,iz)==0.) then
                continue
              else
                scansum(ix,iy,iz) = &
                  x1*scansum(ix,iysh,izsh)+x2*scan(ix,iy,iz)
              endif
            enddo   
          enddo   
          do iy = iy3,iy4,idy
            do ix = 1,nx
              scansum(ix,iy,iz) = scan(ix,iy,iz)
            enddo   
          enddo   
        enddo   
        do iz = iz3,iz4,idz
          do iy = 1,ny
            do ix = 1,nx
              scansum(ix,iy,iz) = scan(ix,iy,iz)
            enddo   
          enddo   
        enddo   
        do iz = nzscan+1,nz
          do iy = iy1,iy4,idy
            do ix = 1,nx
              iysh = iy+ishy
              if((iysh>=1).and.(iysh<=ny) &
                .and.(scansum(ix,iysh,iz)/=0.)) then
                if(scan(ix,iy,iz)/=0.) &
                  scansum(ix,iy,iz) = &
                    x1*scansum(ix,iysh,iz)+x2*scan(ix,iy,iz)
              else
                scansum(ix,iy,iz) = scan(ix,iy,iz)
              endif
            enddo   
          enddo   
        enddo   
        if(nz<mp) then
          nz = nz+1
          iz = nz
          do iy = iy1,iy4,idy
            do ix = 1,nx
              iysh = iy+ishy
              if((iysh>=1).and.(iysh<=ny) &
                .and.(scansum(ix,iysh,iz)/=0.)) then
                if(std(ix,iy)/=0.) &
                  scansum(ix,iy,iz) = &
                    x1*scansum(ix,iysh,iz)+x2*std(ix,iy)
              else
                scansum(ix,iy,iz) = std(ix,iy)
              endif
              scan(ix,iy,iz) = std(ix,iy)
            enddo   
          enddo   
        endif

      endif

        if(im>1) then
          print '(" Plotting summed scan map")'
          call scanmap(scansum,std,nzscan)
        endif

        sumtime = sumtime+beamtime
        ierr = 0
        return
      end


      subroutine subscansky(scan,std,nz,ierr)
!	subtract interpolated sky from scan

        use ius
        use dims
        use modes

        real scan(mx,my,1),std(mx,my),sky(mx,my),dsky(mx,my)
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical savep
        logical logstr
        character(80) yn
        character(40) line

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)

!	if(intsky(1)<0) call subscancorr(scan,std,nz,ierr)
        if(intsky(1)<0) call oldscancorr(scan,std,nz,ierr)

        do iy = 1,ny
          do ix = 1,nx
            sky(ix,iy) = 0.
            dsky(ix,iy) = 0.
          enddo   
        enddo   
        nzsky = 0
        zm = 0.
        dnz = 0.

  104   iz1 = abs(intsky(1))
        iz2 = abs(intsky(2))
        iz3 = intsky(3)
        iz4 = intsky(4)
        if(intsky(1)<0) then
          iz3 = 0
          iz4 = 0
        endif

        if((iz1>=1).and.(iz2>iz1).and.(iz2<=nz)) then

          nzsky = iz2-iz1+1
          zm = (iz1+iz2)/2.
          if((iz3>iz2).and.(iz4>=iz3).and.(iz4<=nz)) then
            nzsky = nzsky+iz4-iz3+1
            zm = ((iz2-iz1+1.)*(iz2+iz1)+(iz4-iz3+1.)*(iz4+iz3)) &
                /(2.*nzsky)
  
            do iz = iz3,iz4
              zsky = iz-zm
              dnz = dnz+zsky**2
              do iy = 1,ny
                do ix = 1,nx
                  if(illum(ix,iy)<=0) cycle
                  sky(ix,iy) = sky(ix,iy)+scan(ix,iy,iz)
                  dsky(ix,iy) = dsky(ix,iy)+zsky*scan(ix,iy,iz)
                enddo   
              enddo   
            enddo   
          endif

          do iz = iz1,iz2
            zsky = iz-zm
            if(dnz>0.) dnz = dnz+zsky**2
            do iy = 1,ny
              do ix = 1,nx
                if(illum(ix,iy)<=0) cycle
                sky(ix,iy) = sky(ix,iy)+scan(ix,iy,iz)
                dsky(ix,iy) = dsky(ix,iy)+zsky*scan(ix,iy,iz)
              enddo   
            enddo   
          enddo   

        elseif(iz1==0) then

!	make postage stamp images for user to find sky

          print '(" Plotting scan images.  Find regions on sky")'
          savep = pause
          pause = .true.
          call scanmap(scan,std,nz)

          print '(" Enter skyint (y) or sky frame numbers (n)? ",$)'
          read '(a80)', yn
          if(logstr(yn,ierr)) then
  118       print '(" Enter skyint(1,2[,3,4]) ",$)'
            read '(a40)', line
            read(unit=line,fmt=*,err=118,end=118) (intsky(i),i=1,2)
            intsky(3) = 0
            intsky(4) = 0
            read(unit=line,fmt=*,err=119,end=119) (intsky(i),i=1,4)
  119       continue
            pause = savep
            goto 104
          else
            print '(" Enter iz on sky, or 0 to finish ",$)'
  110       read *, izsky
            if(izsky<=0) goto 120
            nzsky = nzsky+1
            do iy = 1,ny
              do ix = 1,nx
                if(illum(ix,iy)<=0) cycle
                sky(ix,iy) = sky(ix,iy)+scan(ix,iy,izsky)
              enddo   
            enddo   
            goto 110
          endif
  120     pause = savep

        endif

!	subtract interpolated sky 

        if(nzsky>0) then

          xnzsky = nzsky
          do iy = 1,ny
            do ix = 1,nx
              sky(ix,iy) = sky(ix,iy)/xnzsky
              if(dnz>0.) dsky(ix,iy) = dsky(ix,iy)/dnz
            enddo   
          enddo   

          do iz = 1,nz
            zsky = iz-zm
            do iy = 1,ny
              do ix = 1,nx
                if(dnz>0.) then
                  scan(ix,iy,iz) = scan(ix,iy,iz)-sky(ix,iy) &
                    -zsky*dsky(ix,iy)
                else
                  scan(ix,iy,iz) = scan(ix,iy,iz)-sky(ix,iy)
                endif
              enddo   
            enddo   
          enddo   

          if(nz<mp) then
            iz = nz+1
            do iy = 1,ny
              do ix = 1,nx
                scan(ix,iy,iz) = scan(ix,iy,iz)+sky(ix,iy)
              enddo   
            enddo   
          endif

        endif

        return
      end


      subroutine oldscancorr(scan,std,nz,ierr)
!	subtract correlated sky fluctuations from scan
!	called if intsky(1) < 0

        use ius
        use dims
        use modes

        real scan(mx,my,1),std(mx,my)
        real sky(mx,my),dsky(mp),avg(mx,my)
        real krot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /iwin/ iwin1,iwin2

        print '(" Removing correlated sky noise from scan")'
        write(iupl,'(" Removing correlated sky noise from scan")')

        print '(" Plotting noisy scan map")'
        call scanmap(scan,std,nz)

        iorder1 = abs(intsky(3))
        iorder2 = abs(intsky(4))
        if((iorder1<1).or.(iorder2>norder) &
          .or.(iorder1>iorder2)) then
          iorder1 = 1
          iorder2 = norder
        endif
        it1 = intshift(1)
        it2 = intshift(2)
        if((it1<1).or.(it2>nt).or.(it1>it2)) then
          it1 = 1
          it2 = nt
        endif

!	first guess for sky noise: flat or spatial average of std

        do iy = 1,ny
          do iorder = 1,norder
            if(crossdisp) then
              sum = 0.
              nti = 0
!	only use skyint orders to determine fluctuations in sky
              if((iorder>=iorder1).and.(iorder<=iorder2)) then
                do it = it1,it2
                  ix = (iorder-1)*nt+it
                  if((ix>=1).and.(ix<=nx) &
                    .and.(illum(ix,iy)>0)) then
                    sum = sum+std(ix,iy)
                    nti = nti+1
                  endif
                enddo   
              endif
              if(nti>0) then
                sum = sum/nti
              else
                sum = 0.
              endif
            else
              sum = 1.
            endif
            do it = 1,nt
              ix = (iorder-1)*nt+it
              if((ix>=1).and.(ix<=nx)) then
                if(illum(ix,iy)>0) then
                  sky(ix,iy) = sum
                else
                  sky(ix,iy) = 0.
                endif
              endif
            enddo   
          enddo   
        enddo   

        do iy = 1,ny
          do ix = 1,nx
            if(illum(ix,iy)>0) then
              sum = 0.
              do iz = 1,nz
                sum = sum+scan(ix,iy,iz)
              enddo   
              avg(ix,iy) = sum/nz
            else
              avg(ix,iy) = 0.
            endif
          enddo   
        enddo   

        sumsky = 0.
        do iy = 1,ny
          do iorder = iorder1,iorder2
            do it = it1,it2
              ix = (iorder-1)*nt+it
              if((ix>=1).and.(ix<=nx) &
                .and.(illum(ix,iy)>0)) then
                sumsky = sumsky+sky(ix,iy)**2
              endif
            enddo   
          enddo   
        enddo   

!	find amount of sky in each frame

        sumdsky = 0.
        do iz = 1,nz
          sum = 0.
          do iy = 1,ny
            do iorder = iorder1,iorder2
              do it = it1,it2
                ix = (iorder-1)*nt+it
                if((ix>=1).and.(ix<=nx) &
                  .and.(illum(ix,iy)>0)) then
                  sum = sum+sky(ix,iy)*(scan(ix,iy,iz)-avg(ix,iy))
                endif
              enddo   
            enddo   
          enddo   
          dsky(iz) = sum/sumsky
          sumdsky = sumdsky+dsky(iz)**2
        enddo   

        if(.not.crossdisp) goto 290
!	second guess for sky: correlation of scan frames with dsky

        sumsky = 0.
        do iy = 1,ny
          do iorder = 1,norder
            sum = 0.
            nti = 0
            do it = it1,it2
              ix = (iorder-1)*nt+it
              if((ix>=1).and.(ix<=nx) &
                .and.(illum(ix,iy)>0)) then
                avgi = avg(ix,iy)
                do iz = 1,nz
                  sum = sum+dsky(iz)*(scan(ix,iy,iz)-avgi)
                enddo   
                nti = nti+1
              endif
            enddo   
            if(nti>0) then
              sum = sum/(nti*sumdsky)
            else
              sum = 0.
            endif
            do it = 1,nt
              ix = (iorder-1)*nt+it
              if((ix>=1).and.(ix<=nx)) then
                if(illum(ix,iy)>0) then
                  sky(ix,iy) = sum
                  if((iorder>=iorder1).and.(iorder<=iorder2) &
                    .and.(it>=it1).and.(it<=it2)) &
                    sumsky = sumsky+sum**2
                else
                  sky(ix,iy) = 0.
                endif
              endif
            enddo   
          enddo   
        enddo   

        if(verbose) then
          print '(" Plotting sky noise array")'
          call grey(mx,nx,ny,sky,0.,0.,iwin1)
          call waitasec(pause,.true.,sky,mx,nx,my)
        endif

!	find amount of sky in each frame

        do iz = 1,nz
          sum = 0.
          do iy=1,ny
            do iorder = iorder1,iorder2
              do it = it1,it2
                ix = (iorder-1)*nt+it
                if((ix>=1).and.(ix<=nx).and.(illum(ix,iy)>0)) &
                  sum = sum+sky(ix,iy)*(scan(ix,iy,iz)-avg(ix,iy))
              enddo   
            enddo   
          enddo   
          dsky(iz) = sum/sumsky
        enddo   

  290   print '(" Sky fluctuations:")'
        print '(8f9.5)', (dsky(iz),iz=1,nz)

!	subtract sky from scan frames

        do iz = 1,nz
          do iy = 1,ny
            do ix = 1,nx
              scan(ix,iy,iz) = scan(ix,iy,iz)-dsky(iz)*sky(ix,iy)
            enddo   
          enddo   
        enddo   

        print '(" Plotting cleaned scan map")'

        return
      end


      subroutine subscancorr(scan,std,nz,ierr)
!	subtract correlated sky fluctuations from scan
!	called if intsky(1) < 0
!	I think this version worked for the GC where there was always
!	blank sky along the slit.
!	But it looks like it's hardwired for xd data.

        use ius
        use dims
        use modes

        real scan(mx,my,1),std(mx,my)
        real sky(mx,my),dsky(mp),avg(mx,my)
        integer isky(mp)
        real krot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /iwin/ iwin1,iwin2

        if(norder<=1) then
          print '("*** subscancorr is only implemented", &
                  " for cross-dispersed data")'
          ierr = 9
          return
        endif
        print '(" Removing correlated sky noise from scan")'
        write(iupl,'(" Removing correlated sky noise from scan")')

        print '(" Plotting noisy scan map")'
        call scanmap(scan,std,nz)

        lz = nz
!	the next line gets sky noise from only the first 6 frames
!	nz = 6

        iorder1 = abs(intsky(3))
        iorder2 = abs(intsky(4))
        if((iorder1<1).or.(iorder2<iorder1) &
          .or.(iorder2>norder)) then
          iorder1 = 1
          iorder2 = norder
        endif
        if((intshift(1)>1).and.(intshift(2)>intshift(1)) &
          .and.(intshift(2)<nt)) then
          it1 = intshift(1)
          it2 = intshift(2)
        else
          it1 = 1
          it2 = nt
        endif

!	find amount of sky in each frame from the minimum along the slit
!	of the spectral average

        xnz = nz
        do iy = 1,ny
          do ix = 1,nx
            if(illum(ix,iy)>0) then
              sum = 0.
              do iz = 1,nz
                sum = sum+scan(ix,iy,iz)
              enddo   
              avg(ix,iy) = sum/xnz
            else
              avg(ix,iy) = 0.
            endif
          enddo   
        enddo   

        sum1 = 0.
        do iz = 1,nz
          spatmin = 0.
          do it = it1,it2
            sum = 0.
            nill = 0
            do iorder = iorder1,iorder2
              ix = (iorder-1)*nt+it
              if(illx(ix)<=0) cycle
              do iy = 1,ny
                if(illum(ix,iy)>0) then
                  sum = sum+scan(ix,iy,iz)
                  nill = nill+1
                endif
              enddo   
            enddo   
            if(nill>0) then
              spat = sum/nill
              if((spatmin==0.).or.(spat<spatmin)) then
                spatmin = spat
                iskyi = it
              endif
            endif
          enddo   
          dsky(iz) = spatmin
          sum1 = sum1+spatmin
          isky(iz) = iskyi
        enddo   

        sum1 = sum1/nz
        sum2 = 0.
        do iz = 1,nz
          dsky(iz) = dsky(iz)-sum1
          sum2 = sum2+dsky(iz)**2
        enddo   

        print '(" Sky fluctuations:")'
        print '(8f9.3)', (dsky(iz),iz=1,nz)
        print '(16i4)', (isky(iz),iz=1,nz)

!	sky(ix,iy) = correlation of scan frames with dsky

        do iy = 1,ny
          do iorder = 1,norder
            if(crossdisp) then
              sum = 0.
              nti = 0
              do it = it1,it2
                ix = (iorder-1)*nt+it
                if((ix>=1).and.(ix<=nx) &
                  .and.(illum(ix,iy)>0)) then
                  avgi = avg(ix,iy)
                  do iz = 1,nz
                    sum = sum+dsky(iz)*(scan(ix,iy,iz)-avgi)
                  enddo   
                  nti = nti+1
                endif
              enddo   
              if(nti>0) then
                sum = sum/(nti*sum2)
              else
                sum = 0.
              endif
            else
              sum = 1.
            endif
            do it = 1,nt
              ix = (iorder-1)*nt+it
              if((ix>=1).and.(ix<=nx)) then
                if(illum(ix,iy)>0) then
                  sky(ix,iy) = sum
                else
                  sky(ix,iy) = 0.
                endif
              endif
            enddo   
          enddo   
        enddo   

        if(verbose.and.crossdisp) then
          print '(" Plotting sky noise array")'
          call grey(mx,nx,ny,sky,0.,0.,iwin1)
          call waitasec(pause,.true.,sky,mx,nx,my)
        endif

!	subtract sky from scan frames and save sky

        do iz = 1,nz
          dskyz = dsky(iz)
          do iy = 1,ny
            do ix = 1,nx
              scan(ix,iy,iz) = scan(ix,iy,iz)-dskyz*sky(ix,iy)
            enddo   
          enddo   
        enddo   
        nz = lz
        if(nz>=mp) return
        do iy = 1,ny
          do ix = 1,nx
            scan(ix,iy,nz+1) = sky(ix,iy)
          enddo   
        enddo   

        return
      end


      subroutine scanmap(scan,std,nz)
!	make postage stamp images for scans

        use dims
        use modes

        real scan(mx,my,1),std(mx,my)
        real map(mx,my)
        real krot
        logical good
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /iwin/ iwin1,iwin2
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf

!	collapse spectrally

      if(crossdisp.and.(modetort/=MNONE)) then

        io1 = intshift(3)
        io2 = intshift(4)
        if((io1<1).or.(io2<io1).or.(io2>norder)) then
          io1 = 1
          io2 = norder
        endif
        iy1 = intshift(5)
        iy2 = intshift(6)
        if((iy1<1).or.(iy2<iy1).or.(iy2>ny)) then
          iy1 = 1
          iy2 = ny
        endif
        nymin = (iy2-iy1+2)/2
        do iz = 1,nz
          do ix = 1,nx
            sum = 0.
            swtsum = 0.
            illsum = 0
            do iy = iy1,iy2
              if((illum(ix,iy)==1).and.(std(ix,iy)>0.)) then
                swt = 1./std(ix,iy)**2
                sum = sum+swt*scan(ix,iy,iz)
                swtsum = swtsum+swt
                illsum = illsum+1
              endif
            enddo   
            if(illsum>nymin) then
              avg = sum/swtsum
            else
              avg = 0.
            endif
            if(nz<=my/2) then
              map(ix,2*iz-1) = avg
              map(ix,2*iz) = avg
            else
              map(ix,iz) = avg
            endif
          enddo   
        enddo   
        if(nz<=my/2) then
          iz2 = 2*nz+1
        else
          iz2 = nz+1
        endif
        do iz = 1,iz2
          do it = 1,nt
            ixsum = nt*norder+it
            if(ixsum>nx) cycle
            summap = 0.
            do iorder = io1,io2
              ix = it+(iorder-1)*nt
              summap = summap+map(ix,iz)
            enddo   
            avg = summap/(io2-io1+1)
            map(ixsum,iz) = avg
          enddo   
        enddo   
        if(nz<=my/2) then
          iz1 = 2*nz+1
        else
          iz1 = nz+1
        endif
        do iz = iz1,my
          do ix = 1,mx
            map(ix,iz) = 0.
          enddo   
        enddo   
        if(ask.and.verbose) then
          call grey(mx,mx,my,map,-1.,-1.,iwin1)
        else
          call grey(mx,mx,my,map,0.,0.,iwin1)
        endif
        call waitasec(pause,.true.,map,mx,mx,my)

      elseif((modeinst==MMED).or.(modeinst==MLOW)) then

        ix1 = intshift(5)
        ix2 = intshift(6)
        if((ix1<1).or.(ix2<ix1).or.(ix2>nx)) then
          ix1 = 1
          ix2 = nx
        endif
        nxmin = (ix2-ix1+2)/2
        if((ix1<1).or.(ix2<ix1).or.(ix2>nx)) then
          ix1 = 1
          ix2 = nx
        endif
        do iz = 1,nz
          sumavg = 0.
          illavg = 0
          ip = 0
          do iy = 1,ny
            sum = 0.
            swtsum = 0.
            do ix = ix1,ix2
              if((illum(ix,iy)==1).and.(std(ix,iy)>0.)) then
                swt = 1./std(ix,iy)**2
                sum = sum+swt*scan(ix,iy,iz)
                swtsum = swtsum+swt
              endif
            enddo   
            if(swtsum>0.) then
              avg = sum/swtsum
              sumavg = sumavg+avg
              illavg = illavg+1
            else
              avg = 0.
            endif
            if(nz<=my/2) then
              map(iy,2*iz-1) = avg
              map(iy,2*iz) = avg
            else
              map(iy,iz) = avg
            endif
          enddo   
          if(dosubsky) then
            avg = sumavg/illavg
            do iy = 1,ny
              if(nz<=my/2) then
                if(map(iy,2*iz)/=0.) then
                  map(iy,2*iz-1) = map(iy,2*iz-1)-avg
                  map(iy,2*iz) = map(iy,2*iz)-avg
                endif
              else
                if(map(iy,iz)/=0.) map(iy,iz) = map(iy,iz)-avg
              endif
            enddo   
          endif
        enddo   
        if(nz<=my/2) then
          iz1 = 2*nz+1
        else
          iz1 = nz+1
        endif
        do iz = iz1,my
          do iy = 1,my
            map(iy,iz) = 0.
          enddo   
        enddo   
        if(ask.and.verbose) then
          call grey(my,my,my,map,-1.,-1.,iwin1)
        else
          call grey(my,my,my,map,0.,0.,iwin1)
        endif
        call waitasec(pause,.true.,map,my,my,my)

      else
!	scanning in camera mode?

        do iz = 1,nz
          sumavg = 0.
          illavg = 0
          ip = 0
          do ix = 1,nx
            sum = 0.
            swtsum = 0.
            do iy = 1,ny
              if((illum(ix,iy)==1).and.(std(ix,iy)>0.)) then
                swt = 1./std(ix,iy)**2
                sum = sum+swt*scan(ix,iy,iz)
                swtsum = swtsum+swt
              endif
            enddo   
            if(swtsum>0.) then
              avg = sum/swtsum
              sumavg = sumavg+avg
              illavg = illavg+1
            else
              avg = 0.
            endif
            map(ix,iz) = avg
          enddo   
          if(dosubsky) then
            avg = sumavg/illavg
            do ix = 1,nx
              if(map(ix,iz)/=0.) map(ix,iz) = map(ix,iz)-avg
            enddo   
          endif
        enddo   
        do iz = nz+1,my
          do ix = 1,mx
            map(ix,iz) = 0.
          enddo   
        enddo   
        if(ask.and.verbose) then
          call grey(mx,nx,nz,map,-1.,-1.,iwin1)
        else
          call grey(mx,nx,nz,map,0.,0.,iwin1)
        endif
        call waitasec(pause,.true.,map,mx,nx,nz)

      endif

        ierr = 0
        return
      end


      subroutine scanmap2(scan,scansum,std,nz0)
!	make postage stamp images of scan and scansum

        use dims
        use modes

        real scan(mx,my,1),scansum(mx,my,1),std(mx,my)
        real map(mx,my)
        real krot
        logical good,click
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /iwin/ iwin1,iwin2
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf

        click = ask.and.doshift
        nz2 = min(2*nz0,my)
        nz = nz2/2

!	collapse spectrally

      if(crossdisp.and.(modetort/=MNONE)) then

        io1 = intshift(3)
        io2 = intshift(4)
        if((io1<1).or.(io2<io1).or.(io2>norder)) then
          io1 = 1
          io2 = norder
        endif
        iy1 = intshift(5)
        iy2 = intshift(6)
        if((iy1<1).or.(iy2<iy1).or.(iy2>ny)) then
          iy1 = 1
          iy2 = ny
        endif
        nymin = (iy2-iy1+2)/2
        do iz = 1,nz
          do ix = 1,nx
            sum = 0.
            sums = 0.
            swtsum = 0.
            illsum = 0
            do iy = iy1,iy2
              if((illum(ix,iy)==1).and.(std(ix,iy)>0.)) then
                swt = 1./std(ix,iy)**2
                sum = sum+swt*scan(ix,iy,iz)
                sums = sums+swt*scansum(ix,iy,iz)
                swtsum = swtsum+swt
                illsum = illsum+1
              endif
            enddo   
            if(illsum>nymin) then
              avg = sum/swtsum
              avgs = sums/swtsum
            else
              avg = 0.
              avgs = 0.
            endif
            map(ix,iz) = avg
            if((nz+iz)<=my) map(ix,nz+iz) = avgs
          enddo   
        enddo   
        do iz = 1,nz2
          do it = 1,nt
            ixsum = nt*norder+it
            if(ixsum>nx) cycle
            summap = 0.
            do iorder = io1,io2
              ix = it+(iorder-1)*nt
              summap = summap+map(ix,iz)
            enddo   
            map(ixsum,iz) = summap/(io2-io1+1)
          enddo   
        enddo   
        do iz = nz2+1,my
          do ix = 1,mx
            map(ix,iz) = 0.
          enddo   
        enddo   
        if(ask.and.verbose) then
          call grey(mx,nx,ny,map,-1.,-1.,iwin1)
        else
          call grey(mx,nx,ny,map,0.,0.,iwin1)
        endif
        call waitasec(click,.true.,map,mx,nx,ny)

      elseif((modeinst==MMED).or.(modeinst==MLOW)) then

        ix1 = intshift(5)
        ix2 = intshift(6)
        if((ix1<1).or.(ix2<ix1).or.(ix2>nx)) then
          ix1 = 1
          ix2 = nx
        endif
        do iz = 1,nz
          sumavg = 0.
          sumavgs = 0.
          illavg = 0
          ip = 0
          do iy = 1,ny
            sum = 0.
            sums = 0.
            swtsum = 0.
            do ix = ix1,ix2
              if((illum(ix,iy)==1).and.(std(ix,iy)>0.)) then
                swt = 1./std(ix,iy)**2
                sum = sum+swt*scan(ix,iy,iz)
                sums = sums+swt*scansum(ix,iy,iz)
                swtsum = swtsum+swt
              endif
            enddo   
            if(swtsum>0.) then
              avg = sum/swtsum
              avgs = sums/swtsum
              sumavg = sumavg+avg
              sumavgs = sumavgs+avgs
              illavg = illavg+1
            else
              avg = 0.
              avgs = 0.
            endif
            map(iy,iz) = avg
            if((nz+iz)<=my) map(iy,nz+iz) = avgs
          enddo   
          if(dosubsky) then
            avg = sumavg/illavg
            avgs = sumavgs/illavg
            do iy = 1,ny
              if(map(iy,iz)/=0.) map(iy,iz) = map(iy,iz)-avg
              if(((nz+iz)<=my).and.(map(iy,nz+iz)/=0.)) &
                map(iy,nz+iz) = map(iy,nz+iz)-avgs
            enddo   
          endif
        enddo   
        do iz = nz2+1,my
          do iy = 1,my
            map(iy,iz) = 0.
          enddo   
        enddo   
        if(ask.and.verbose) then
          call grey(mx,nx,ny,map,-1.,-1.,iwin1)
        else
          call grey(mx,nx,ny,map,0.,0.,iwin1)
        endif
        call waitasec(click,.true.,map,mx,nx,ny)

      else

        do iz = 1,nz
          sumavg = 0.
          sumavgs = 0.
          illavg = 0
          ip = 0
          do ix = 1,nx
            sum = 0.
            sums = 0.
            swtsum = 0.
            do iy = 1,ny
              if((illum(ix,iy)==1).and.(std(ix,iy)>0.)) then
                swt = 1./std(ix,iy)**2
                sum = sum+swt*scan(ix,iy,iz)
                sums = sums+swt*scansum(ix,iy,iz)
                swtsum = swtsum+swt
              endif
            enddo   
            if(swtsum>0.) then
              avg = sum/swtsum
              avgs = sums/swtsum
              sumavg = sumavg+avg
              sumavgs = sumavgs+avgs
              illavg = illavg+1
            else
              avg = 0.
              avgs = 0.
            endif
            map(ix,iz) = avg
            if((nz+iz)<=my) map(ix,nz+iz) = avgs
          enddo   
          if(dosubsky) then
            avg = sumavg/illavg
            avgs = sumavgs/illavg
            do ix = 1,nx
              if(map(ix,iz)/=0.) map(ix,iz) = map(ix,iz)-avg
              if(((nz+iz)<=my).and.(map(ix,nz+iz)/=0.)) &
                map(ix,nz+iz) = map(ix,nz+iz)-avgs
            enddo   
          endif
        enddo   
        do iz = nz2+1,my
          do ix = 1,mx
            map(ix,iz) = 0.
          enddo   
        enddo   
        if(ask.and.verbose) then
          call grey(mx,nx,ny,map,-1.,-1.,iwin1)
        else
          call grey(mx,nx,ny,map,0.,0.,iwin1)
        endif
        call waitasec(click,.true.,map,mx,nx,ny)

      endif

        ierr = 0
        return
      end


      subroutine writeparms(iunit,junit,ierr)
!	write reduction parameters to reduced header

        real krot,lskrot
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        character(60) comment,NOC

        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /parms/ thrfac,spikefac,stdfac,satval,darkval,xnlin, &
                ynlin,znlin,cloud,bounce,dtcal,dtsky
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /tort0/ hrfl0,xdfl0,xdmrdgr,xdlrdgr,xdmr0,xdlr0,fred0, &
                xdkrot,lskrot

        write(iunit,'("cardmode= ",i1)',err=90) modecard
        write(iunit,'("thrfac  = ",f8.3)',err=90) thrfac
        write(iunit,'("spikefac= ",f8.3)',err=90) spikefac
        write(iunit,'("stdfac  = ",f8.3)',err=90) stdfac
        write(iunit,'("satval  = ",f8.0)',err=90) satval
        write(iunit,'("xnlin   = ",f8.0)',err=90) xnlin
        write(iunit,'("ynlin   = ",f8.0)',err=90) ynlin
        write(iunit,'("znlin   = ",f8.0)',err=90) znlin
        write(iunit,'("cloud   = ",f8.3)',err=90) cloud
        write(iunit,'("bounce  = ",f8.3)',err=90) bounce
        write(iunit,'("dtcal   = ",f8.3)',err=90) dtcal
        write(iunit,'("dtsky   = ",f8.3)',err=90) dtsky
        write(iunit,'("slitrot = ",f10.6)',err=90) slitrot
        write(iunit,'("krot    = ",f10.6)',err=90) krot
        write(iunit,'("detrot  = ",f10.6)',err=90) detrot
        write(iunit,'("hrfl    = ",f8.3)',err=90) hrfl
        write(iunit,'("hrr     = ",f10.6)',err=90) abs(hrr)
        write(iunit,'("hrg     = ",f10.6)',err=90) hrg
        write(iunit,'("hrdgr   = ",f10.6)',err=90) hrdgr
        write(iunit,'("xdfl    = ",f10.3)',err=90) xdfl
        write(iunit,'("xdr     = ",f8.3)',err=90) xdr
        write(iunit,'("xdg     = ",f10.6)',err=90) xdg
        write(iunit,'("xddgr   = ",f10.6)',err=90) xddgr
        write(iunit,'("brl     = ",f10.6)',err=90) brl
        write(iunit,'("x0brl   = ",f10.6)',err=90) x0brl
        write(iunit,'("y0brl   = ",f10.6)',err=90) y0brl
        write(iunit,'("fred    = ",f10.6)',err=90) fred0

        if(dofits) then
          NOC = ' '
          comment = "The following 24 lines are pipe reduction parameters"
          call fithcomm('COMMENT ',comment,junit)
          call fithint('CARDMODE',modecard,NOC,junit)
          call fithreal('THRFAC  ',thrfac,NOC,junit)
          call fithreal('SPIKEFAC',spikefac,NOC,junit)
          call fithreal('STDFAC  ',stdfac,NOC,junit)
          call fithreal('SATVAL  ',satval,NOC,junit)
          call fithreal('XNLIN   ',xnlin,NOC,junit)
          call fithreal('CLOUD   ',cloud,NOC,junit)
          comment = "grating bounce was removed if != 0"
          call fithreal('BOUNCE  ',bounce,comment,junit)
          call fithreal('SLITROT ',slitrot,NOC,junit)
          call fithreal('KROT    ',krot,NOC,junit)
          call fithreal('DETROT  ',detrot,NOC,junit)
          call fithreal('HRFL    ',hrfl,NOC,junit)
          comment = "may have changed after header write"
          call fithreal('HRR     ',abs(hrr),comment,junit)
          call fithreal('HRG     ',hrg,NOC,junit)
          call fithreal('HRDGR   ',hrdgr,NOC,junit)
          call fithreal('XDFL    ',xdfl,NOC,junit)
          call fithreal('XDR     ',xdr,NOC,junit)
          call fithreal('XDG     ',xdg,NOC,junit)
          call fithreal('XDDGR   ',xddgr,NOC,junit)
          call fithreal('BRL     ',brl,NOC,junit)
          call fithreal('X0BRL   ',x0brl,NOC,junit)
          call fithreal('Y0BRL   ',y0brl,NOC,junit)
          call fithreal('FRED    ',fred0,NOC,junit)
        endif

        ierr = 0
        return

   90   ierr = 8
        print '(" Error writing parms to header")'
        return
      end



      subroutine storearr(arr,std,outfile,fitsfile,ierr)
!	store reduced camera, long-slit, or untorted data array
!	to both pipe and reduced fits format files
!	putting iy=1-4 in fits extensions

        use ius
        use dims
        use modes

        real arr(mx,my),std(mx,my)
        real temparr(mx)
        character(32) outfile
        character(48) fitsfile
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical doflip,doflop
        character(60) bunits,comment

        common /byteflip/ doflip,doflop
        common /nn/ nx,ny,nc,norder,ns,nt
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf

!	if(modeobs==MCAMERA) then
          nz = 1
!	else
!	  nz = 2
!	endif

        write(iuredh,'("iunits  = erg/s cm2 sr cm-1")')
        write(iuredh,'("funits  = Jy")')
        write(iuredh,'("bitpix  = -32")')
        write(iuredh,'("nx      = ",i6)') nx
        write(iuredh,'("ny      = ",i6)') ny
        write(iuredh,'("nz      = ",i6)') nz

        line = 4*nx
        if(rdfits) then
          open(unit=iuredd,status='scratch',access='direct', &
            recl=line,err=90)
        else
          open(unit=iuredd,file=outfile,access='direct', &
            recl=line,err=90)
        endif
        do iy = 1,ny
          if(doflop) then
            call flipr4(arr(1,iy),temparr,nx)
            write(unit=iuredd,rec=iy) (temparr(ix),ix=1,nx)
          else
            write(unit=iuredd,rec=iy) (arr(ix,iy),ix=1,nx)
          endif
        enddo   
        if(nz==2) then
          do iy = 1,ny
            jy = iy+ny
            if(doflop) then
              call flipr4(std(1,iy),temparr,nx)
              write(unit=iuredd,rec=jy) (temparr(ix),ix=1,nx)
            else
              write(unit=iuredd,rec=jy) (std(ix,iy),ix=1,nx)
            endif
          enddo   
        endif

        close(unit=iuredd)
        close(unit=iuredh)

        ierr = 0
        if(.not.dofits) return
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

        if(nz==2) then
          call redxtend(iufred,nx,ny-4,0,-3)
          do iy = 5,ny
            call fitsrdat(iufred,std(1,iy),nx)
          enddo   
          call fitsrpad(iufred)
        endif

        call fitsclos(iufred)
        return

   90   print '(" Error writing array")'
        ierr = 9
        close(unit=iuredh)
        return

      end


      subroutine storescan(scan,flat,std,nz,outfile,fitsfile, &
        steptime,nssum,ierr)
!	store scan array

        use ius
        use dims
        use modes
        use consts

        real scan(mx,my,1),flat(mx,my),std(mx,my)
        real frame(mx,my),zero(mx),freq(ms,4),plot(ms,6),black(ms)
        real temparr(mx),flfreq(ms,4),flzero(ms),tbl(ms,6)
        real nodpa,lores,kmirror,krot
        character(32) outfile,nomon
        character(48) fitsfile
        character(60) bunits,comment
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical doflip,doflop
        logical iuopen,firstfile

        common /byteflip/ doflip,doflop
        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdri/ nframe,nsum,nwrite,nnod,nscan,nspec,nspat,nbitpix, &
                nsky,lorder
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /extparms/ modeext,intext(4),intsky(4),intshift(6), &
                intill(4)
        common /iwin/ iwin1,iwin2
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime

        if(dosum) then
!	storing sum over multiple files
          iuhd = iusumh
          iufhd = iufish
        else
!	storing one or more scans from one file
          iuhd = iuredh
          iufhd = iufith
        endif
        ierr = 0
        inquire(unit=iuhd,opened=iuopen)
        if(.not.iuopen) then
          if(rdfits) then
            open(unit=iuhd,status='scratch')
          else
            print '("*** storescan called with hdr not opened",i4)', iuhd
            stop
          endif
        endif

!	slitoff must be copied from extract
        slitoff = -0.3

        if(dosum) then
          if(verbose) then
!	plotting scansum
            print '(" Plotting summed scan frames")'
            iwinn = iwin1
            do iz = 1,nz
              do iy = 1,ny
                do ix = 1,nx
                  frame(ix,iy) = scan(ix,iy,iz)
                enddo   
              enddo   
              call grey(mx,nx,ny,frame,0.,0.,iwinn)
              if(pause.and.verbose) then
                call waitasec(pause,.true.,frame,mx,nx,ny)
              else
                call mswait(200,xxx)
              endif
              iwinn = -1
            enddo   
          endif
          print '(" Plotting summed map")'
!	  call scanmap(scan,flat,nz)
          call scanmap(scan,std,nz)
        endif

        firstfile = (sumtime<1.5*nscan*beamtime)
        if((modeinst==MMED).or.(modeinst==MLOW) &
          .or.(modetort==MNONE)) then
!	long-slit or untorted data
          if(firstfile) then
            if(modetort==MNONE) then
              dlnw = 0.
            else
              dlnw = pixelwd/(2.*xdr*xdfl)
            endif
            wnomin = wno0*sexp(dlnw*(1-nx/2))
            wnomax = wno0*sexp(dlnw*(nx/2))
            dvpix = -CLIGHT*dlnw
            rslit = pltscl/(dlnw*slitwid)
            rdiff = 13.5*xdr*wno0
            resolv = 1.0/sqrt(1.0/rslit**2+1.0/rdiff**2)
            write(iuhd,'("wno0    = ",f10.4)') wno0
            write(iuhd,'("wnomin  = ",f10.4)') wnomin
            write(iuhd,'("wnomax  = ",f10.4)') wnomax
            write(iuhd,'("dvpix   = ",f11.6)') dvpix
            write(iuhd,'("resolv  = ",f10.3)') resolv
            write(iuhd,'("sumtime = ",f10.3)') steptime
            write(iuhd,'("iunits  = erg/s cm2 sr cm-1")')
            write(iuhd,'("bitpix  = -32")')
            write(iuhd,'("nx      = ",i6)') nx
            write(iuhd,'("ny      = ",i6)') ny
            write(iuhd,'("nz      = ",i6)') nz
            if(dofits) then
              comment = 'wavenumber at pixel 128 (cm-1)'
              call fithreal('WNO0    ',wno0,comment,iufhd)
              comment = 'minimum wavenumber (cm-1)'
              call fithreal('WNO_MIN ',wnomin,comment,iufhd)
              comment = 'maximum wavenumber (cm-1)'
              call fithreal('WNO_MAX ',wnomax,comment,iufhd)
              comment = 'velocity interval between pixels (km/s)'
              call fithreal('DVPIX   ',dvpix,comment,iufhd)
              comment =  'approximate resolving power'
              call fithreal('RESOLV  ',resolv,comment,iufhd)
              comment = 'summed time per step (s)'
              call fithreal('SUMTIME ',steptime,comment,iufhd)
              bunits = 'erg/(s cm2 sr cm-1)'
              comment = 'cgs intensity units'
              call fithchar('BUNIT   ',bunits,comment,iufhd)
            endif
          else
            write(iuhd,'("sumtime = ",f10.3)') sumtime
            if(dofits) then
              comment = 'summed time per step (s)'
              call fithreal('SUMTIME ',steptime,comment,iufhd)
            endif
          endif
          if(dofits) then
            if(index(fitsfile,' ')>1) then
              call redstart(nx,ny-4,nz,fitsfile,iufhd,ierr)
            else
              ierr = 9
            endif
            if(ierr>0) print '("*** Error opening fits file")'
          endif

          iy1 = max(1,intext(1))
          iy2 = min(nz,intext(2))
          iz1 = max(1,intext(3))
          iz2 = min(ny,intext(4))
          if(iy2<iy1) iy2 = ny
          if(iz2<iz1) iz2 = nz
          sqrtn = sqrt(float(nssum))
          do ix = 1,nx
            x = ix-nx/2
            freq(ix,1) = wno0*sexp(dlnw*x)
            cardsum = 0.
            blksum = 0.
            stdsum = 0.
            isum = 0
            objsum = 0.
            msum = 0
            do iy = 1,ny
              if(cards(ix,iy,1)>0.) then
                cardsum = cardsum+cards(ix,iy,1)
                blksum = blksum+cards(ix,iy,3)
                stdsum = stdsum+std(ix,iy)
                isum = isum+1
              endif
              if((iy>=iy1).and.(iy<=iy2)) then
                do iz = iz1,iz2
                  objsum = objsum+scan(ix,iy,iz)
                  msum = msum+1
                enddo   
              endif
            enddo   
            if(isum>0) then
              freq(ix,3) = stdsum/(sqrtn*isum)
              freq(ix,4) = cardsum/isum
              black(ix) = blksum/isum
            else
              freq(ix,3) = 0.
              freq(ix,4) = 0.
              black(ix) = 0.
            endif
            if(msum>0) then
              freq(ix,2) = objsum/msum
            else
              freq(ix,2) = 0.
            endif
          enddo   
!	wait to make table until after saving data
!          if(dofits.and.(ierr==0)) then
!            do iy = 1,4
!              call redxtend(iufred,nx,1,0,iy)
!              call fitsrdat(iufred,freq(1,iy),nx)
!              call fitsrpad(iufred)
!            enddo   
!          endif

          if(doflop) then
            do it = 1,4
              call flipr4(freq(1,it),flfreq(1,it),nx)
            enddo   
          endif

          line = 4*nx
          if(rdfits) then
            open(unit=iuredd,status='scratch',access='direct', &
              recl=line,err=90)
          else
            open(unit=iuredd,file=outfile,access='direct', &
              recl=line,err=90)
          endif
          irec = 0
          do iz = 1,nz
            do it = 1,4
              irec = irec+1
              if(doflop) then
                write(unit=iuredd,rec=irec) (flfreq(ix,it),ix=1,nx)
              else
                write(unit=iuredd,rec=irec) (freq(ix,it),ix=1,nx)
              endif
            enddo   
!            if(dofits.and.(ierr==0)) &
!              call redxtend(iufred,nx,ny-4,iz,0)
            do iy = 5,ny
              irec = irec+1
              if(doflop) then
                call flipr4(scan(1,iy,iz),temparr,nx)
                write(unit=iuredd,rec=irec) (temparr(ix),ix=1,nx)
              else
                write(unit=iuredd,rec=irec) (scan(ix,iy,iz),ix=1,nx)
              endif
              if(dofits.and.(ierr==0)) &
                call fitsrdat(iufred,scan(1,iy,iz),nx)
            enddo   
          enddo   
          if(dofits.and.(ierr==0)) call fitsrpad(iufred)

          if(dofits.and.(ierr==0)) then
            do ix = 1,nx
              do iy = 1,4
                tbl(ix,iy) = freq(ix,iy)
              enddo
              tbl(ix,5) = 1.0e+04/freq(ix,1)
            enddo
            call tblxtend(iufred,nx,5,tbl)
          endif

          iz = nz
          if(nz<mp) then
!	store sky frame
            iz = nz+1
            do it = 1,4
              irec = irec+1
              write(unit=iuredd,rec=irec) (flfreq(ix,it),ix=1,nx)
            enddo   
            if(dofits.and.(ierr==0)) &
              call redxtend(iufred,nx,ny-4,0,-1)
            do iy = 5,ny
              irec = irec+1
              if(doflop) then
                call flipr4(scan(1,iy,iz),temparr,nx)
                write(unit=iuredd,rec=irec) (temparr(ix),ix=1,nx)
              else
                write(unit=iuredd,rec=irec) (scan(ix,iy,iz),ix=1,nx)
              endif
              if(dofits.and.(ierr==0)) &
                call fitsrdat(iufred,scan(1,iy,iz),nx)
            enddo   
            if(dofits.and.(ierr==0)) call fitsrpad(iufred)
          endif
          if(nz+1<mp) then
!	store noise frame
            iz = nz+2
            do it = 1,4
              irec = irec+1
              write(unit=iuredd,rec=irec) (flfreq(ix,it),ix=1,nx)
            enddo   
            if(dofits.and.(ierr==0)) &
              call redxtend(iufred,nx,ny-4,0,-2)
            do iy = 5,ny
              irec = irec+1
              if(doflop) then
                call flipr4(scan(1,iy,iz),temparr,nx)
                write(unit=iuredd,rec=irec) (temparr(ix),ix=1,nx)
              else
                write(unit=iuredd,rec=irec) (scan(ix,iy,iz),ix=1,nx)
              endif
              if(dofits.and.(ierr==0)) &
                call fitsrdat(iufred,scan(1,iy,iz),nx)
            enddo   
            if(dofits.and.(ierr==0)) call fitsrpad(iufred)
          endif
          if(iz>nz) then
            write(iuhd,'("skyframe= T")')
          else
            write(iuhd,'("skyframe= F")')
          endif

          if(domon) then
            dopp = (1.+radvel/CLIGHT)
            do ix = 1,nx
              plot(ix,1) = dopp*freq(ix,1)
              plot(ix,2) = freq(ix,2)
              plot(ix,3) = freq(ix,3)
              plot(ix,4) = freq(ix,4)
              plot(ix,5) = black(ix)
            enddo   
            nplt = 1
            scale = -1.
            if(ask) scale = 0.
            print '(" Plotting summed spectrum")'
            call plots(ms,nx,nplt,plot,scale,outfile,iwin2)
            call waitasec(pause,.false.,plot,ms,nx,ny)
          endif

        else
!	cross-dispersed data
          if(doflop.and.(ny>mx)) &
            print'("*** This program is about to crash!")'

          do iy = 1,ny
            zero(iy) = 0.
          enddo
          if(doflop) call flipr4(zero,flzero,ny)

          dw = 0.5/(sqrt(hrr**2/(1.+hrr**2))*hrdgr)
!          dlnw = pixelwd/(2.*abs(hrr)*hrfl)
          dlnw = pixelwd/(2.*abs(hrr)*(1.-slitoff/20.)*hrfl)
          if(firstfile) then
            wnoi = wno0+(1-(norder+1)/2)*dw
            wnomin = wnoi*sexp(dlnw*(1-ny/2))
            wnoi = wno0+(norder-(norder+1)/2)*dw
            wnomax = wnoi*sexp(dlnw*(ny/2))
            dvpix = -CLIGHT*dlnw
            rslit = pltscl/(dlnw*slitwid)
            rdiff = 18.0*abs(hrr)*wno0
            resolv = 1.0/sqrt(1.0/rslit**2+1.0/rdiff**2)
            write(iuhd,'("wno0    = ",f10.4)') wno0
            write(iuhd,'("wnomin  = ",f10.4)') wnomin
            write(iuhd,'("wnomax  = ",f10.4)') wnomax
            write(iuhd,'("dvpix   = ",f11.6)') dvpix
            write(iuhd,'("resolv  = ",f10.3)') resolv
            write(iuhd,'("sumtime = ",f10.3)') steptime
            write(iuhd,'("iunits  = erg/s cm2 sr cm-1")')
            write(iuhd,'("bitpix  = -32")')
            write(iuhd,'("nx      = ",i6)') ns
            write(iuhd,'("ny      = ",i6)') nt+4
            write(iuhd,'("nz      = ",i6)') nz
            if(dofits) then
              comment = 'wavenumber at pixel 128 of middle order'
              call fithreal('WNO0    ',wno0,comment,iufhd)
              comment = 'minimum wavenumber (cm-1)'
              call fithreal('WNO_MIN ',wnomin,comment,iufhd)
              comment = 'maximum wavenumber (cm-1)'
              call fithreal('WNO_MAX ',wnomax,comment,iufhd)
              comment =  'delta v per pixel (km/s)'
              call fithreal('DVPIX   ',dvpix,comment,iufhd)
              comment =  'approximate resolving power'
              call fithreal('RESOLV  ',resolv,comment,iufhd)
              comment = 'summed time per step (s)'
              call fithreal('SUMTIME ',steptime,comment,iufhd)
              bunits = 'erg/(s cm2 sr cm-1)'
              comment = 'cgs intensity units'
              call fithchar('BUNIT   ',bunits,comment,iufhd)
            endif
          else
            write(iuhd,'("sumtime = ",f10.3)') sumtime
            if(dofits) then
              comment = 'summed time per step (s)'
              call fithreal('SUMTIME ',steptime,comment,iufhd)
            endif
          endif
          if(dofits) then
            if(index(fitsfile,' ')>1) then
              call redstart(ns,nt,nz,fitsfile,iufhd,ierr)
            else
              ierr = 9
            endif
            if(ierr>0) print '("*** Error opening fits file")'
          endif

          sqrtn = sqrt(float(nssum))
          dw = 0.5/(sqrt(hrr**2/(1.+hrr**2))*hrdgr)
          if((intext(1)>1).and.(intext(2)>intext(1)) &
            .and.(intext(2)<nt)) then
            it1 = intext(1)
            it2 = intext(2)
          else
            it1 = 1
            it2 = nt
          endif
          if((intext(3)>0).and.(intext(4)>=intext(3)) &
            .and.(intext(4)<=nz)) then
            iz1 = intext(3)
            iz2 = intext(4)
          else
            iz1 = 1
            iz2 = nz
          endif

!	calculate first four lines
          do iorder = 1,norder
            wnoi = wno0+(iorder-(norder+1)/2)*dw
            do iy = 1,ny
              is = (iorder-1)*ny+iy
              y = iy-ny/2
              freq(is,1) = wnoi*sexp(dlnw*y)
              freq(is,2) = 0.
              cardsum = 0.
              blksum = 0.
              stdsum = 0.
              isum = 0
              do it = it1,it2
                ix = (iorder-1)*nt+it
                if((ix>=1).and.(ix<=nx).and.(cards(ix,iy,1)>0.)) &
                  then
                  cardsum = cardsum+cards(ix,iy,1)
                  blksum = blksum+cards(ix,iy,3)
                  stdsum = stdsum+std(ix,iy)
                  isum = isum+1
                endif
              enddo   
              if(isum>0) then
                freq(is,3) = stdsum/(sqrtn*isum)
                freq(is,4) = cardsum/isum
                black(is) = blksum/isum
              else
                freq(is,3) = 0.
                freq(is,4) = 0.
                black(is) = 0.
              endif
              dopp = (1.+radvel/CLIGHT)
              plot(is,1) = dopp*freq(is,1)
              plot(is,2) = 0.
              plot(is,3) = freq(is,3)
              plot(is,4) = freq(is,4)
              plot(is,5) = black(is)
            enddo   
          enddo   
!          if(dofits.and.(ierr==0)) then
!            do it = 1,4
!              call redxtend(iufred,ns,1,0,it)
!              call fitsrdat(iufred,freq(1,it),ns)
!              call fitsrpad(iufred)
!            enddo   
!          endif


!	extract fluxes, rearrange scan frames in spec order
!	and write out 4+nt lines
          line = 4*ny
          if(rdfits) then
            open(unit=iuredd,status='scratch',access='direct',recl=line, &
               err=90)
          else
            open(unit=iuredd,file=outfile,access='direct',recl=line, &
               err=90)
          endif
          irec = 0
          do iz = 1,nz
            do iorder = 1,norder
              wnoi = wno0+(iorder-(norder+1)/2)*dw
              do iy = 1,ny
                is = (iorder-1)*ny+iy
                y = iy-ny/2
                sum = 0.
                isum = 0
                do it = it1,it2
                  ix = (iorder-1)*nt+it
                  if((ix>=1).and.(ix<=nx) &
                    .and.(cards(ix,iy,1)>0.)) then
                    sum = sum+scan(ix,iy,iz)
                    isum = isum+1
                  endif
                enddo   
                if(isum>0) then
                  freq(is,2) = sum/isum
                  if((iz>=iz1).and.(iz<=iz2)) &
                    plot(is,2) = plot(is,2)+sum/isum
                else
                  freq(is,2) = 0.
                endif
              enddo   
            enddo   
            do it = 1,4
              do iorder = 1,norder
                is0 = (iorder-1)*ny
                do iy = 1,ny
                  is = is0+iy
                  temparr(iy) = freq(is,it)
                enddo   
                irec = irec+1
                if(doflop) call flipr4(temparr,temparr,ny)
                write(unit=iuredd,rec=irec) (temparr(iy),iy=1,ny)
              enddo   
            enddo   
!            if(dofits.and.(ierr==0)) call redxtend(iufred,ns,nt,iz,0)
            do it = 1,nt
              do iorder = 1,norder
                ix = (iorder-1)*nt+it
                irec = irec+1
                if((ix>=1).and.(ix<=nx)) then
                  do iy = 1,ny
                    temparr(iy) = scan(ix,iy,iz)
                  enddo   
                  if(dofits.and.(ierr==0)) &
                    call fitsrdat(iufred,temparr,ny)
                  if(doflop) call flipr4(temparr,temparr,ny)
                  write(unit=iuredd,rec=irec) (temparr(iy),iy=1,ny)
                else
                  if(dofits.and.(ierr==0)) &
                    call fitsrdat(iufred,zero,ny)
                  write(unit=iuredd,rec=irec) (flzero(iy),iy=1,ny)
                endif
              enddo   
            enddo   
          enddo   
          if(dofits.and.(ierr==0)) call fitsrpad(iufred)

          if(dofits.and.(ierr==0)) then
            do is = 1,ns
              do it = 1,4
                tbl(is,it) = freq(is,it)
              enddo
              tbl(is,5) = 1.0e+04/freq(is,1)
            enddo
            call tblxtend(iufred,ns,5,tbl)
          endif

!	write out sky and noise frames
          iz1 = nz+1
          iz2 = min(nz+2,mp)
          if(iz2>nz) then
            write(iuhd,'("skyframe= T")')
          else
            write(iuhd,'("skyframe= F")')
          endif
          do iz = iz1,iz2
            do it = 1,4
              do iorder = 1,norder
                is0 = (iorder-1)*ny
                do iy=1,ny
                  is = is0+iy
                  temparr(iy) = freq(is,it)
                enddo   
                irec = irec+1
                if(doflop) call flipr4(temparr,temparr,ny)
                write(unit=iuredd,rec=irec) (temparr(iy),iy=1,ny)
              enddo   
            enddo   
            if(dofits.and.(ierr==0)) &
              call redxtend(iufred,ns,nt,0,nz-iz)
            do it = 1,nt
              do iorder = 1,norder
                ix = (iorder-1)*nt+it
                irec = irec+1
                if((ix>=1).and.(ix<=nx)) then
                  do iy = 1,ny
                    temparr(iy) = scan(ix,iy,iz)
                  enddo   
                  if(dofits.and.(ierr==0)) &
                    call fitsrdat(iufred,temparr,ny)
                  if(doflop) call flipr4(temparr,temparr,ny)
                  write(unit=iuredd,rec=irec) (temparr(iy),iy=1,ny)
                else
                  if(dofits.and.(ierr==0)) &
                    call fitsrdat(iufred,zero,ny)
                  write(unit=iuredd,rec=irec) (flzero(iy),iy=1,ny)
                endif
              enddo   
            enddo   
            if(dofits.and.(ierr==0)) call fitsrpad(iufred)
          enddo   

          if(domon) then
            scale = -1.
            if(ask) scale = 0.
            nplt = norder
            print '(" Plotting extracted spectrum")'
            call plots(ms,nx,nplt,plot,scale,outfile,iwin2)
            call waitasec(pause,.false.,plot,ms,ms,nx)
          endif

        endif

        close(unit=iuredd)
        if(dofits.and.(ierr==0)) call fitsclos(iufred)
        if(.not.dosum) close(unit=iuredh)
        ierr = 0
        return

   90   print '(" Error opening output data file")'
        if(iuhd==iuredh) close(unit=iuredh)
        ierr = 9
        return

      end


      subroutine storespec(spec,outfile,fitsfile,ierr)
!	store xd spec array

        use ius
        use dims
        use modes

        real spec(ms,mt),tbl(ms,6)
        real temparr(ms)
        real krot
        character(32) outfile
        character(48) fitsfile
        character(60) comment,bunits
        logical baddata,hdopen
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical doflip,doflop

        common /byteflip/ doflip,doflop
        common /nn/ nx,ny,nc,norder,ns,nt
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl

        if((modeinst==MHIMED).or.(modeinst==MHILOW)) then
          ns = ny*norder
          if(ns>ms) ns = ms
          lt = nt+4
          if(lt>mt) then
            lt = mt
            print '("***Why was nt > mt?")'
          endif
        else
          print '("***Why was storespec called?")'
          ns = nx
          lt = ny
        endif
        ierr = 0
        write(iuredh,'("iunits  = erg/s cm2 sr cm-1")')
        write(iuredh,'("funits  = Jy")')
        write(iuredh,'("bitpix  = -32")')
        write(iuredh,'("nx      = ",i6)') ns
        write(iuredh,'("ny      = ",i6)') lt
        write(iuredh,'("nz      = 1")')

        line = 4*ns
        if(rdfits) then
          open(unit=iuredd,status='scratch',access='direct', &
            recl=line,err=90)
        else
          open(unit=iuredd,file=outfile,access='direct',recl=line,err=90)
        endif
        do it = 1,lt
          if(doflop) then
            call flipr4(spec(1,it),temparr,ns)
            write(unit=iuredd,rec=it) (temparr(is),is=1,ns)
          else
            write(unit=iuredd,rec=it) (spec(is,it),is=1,ns)
          endif
        enddo   

        close(unit=iuredd)
        close(unit=iuredh)

        ierr = 0
        if(.not.dofits) return
        bunits = 'erg/s cm2 sr cm-1'
        comment = 'cgs intensity units'
        call fithchar('BUNIT   ',bunits,comment,iufith)
        if(index(fitsfile,' ')==1) return
        call redstart(ns,nt,1,fitsfile,iufith,ierr)
        if(ierr/=0) then
          print '("*** redstart returned",i4)', ierr
          return
        endif

        do it = 5,lt
          call fitsrdat(iufred,spec(1,it),ns)
        enddo   
        call fitsrpad(iufred)

!        do it = 1,4
!          call redxtend(iufred,ns,1,0,it)
!          call fitsrdat(iufred,spec(1,it),ns)
!          call fitsrpad(iufred)
!        enddo   

        do is = 1,ns
          do it = 1,4
            tbl(is,it) = spec(is,it)
          enddo
          tbl(is,5) = 1.0e+04/spec(is,1)
        enddo
        call tblxtend(iufred,ns,5,tbl)

        call fitsclos(iufred)

        return

   90   print '(" Error writing spec")'
        ierr = 9
        close(unit=iuredh)
        return

      end


      subroutine storesum(outfile,fitsfile,ierr)
!	store summed long-slit or xd spectrum

        use ius
        use dims
        use modes
        use consts

        character(32) outfile
        character(48) fitsfile
        character(60) comment,bunits
        real nodpa,krot,lores,kmirror
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        logical doflip,doflop
        real temparr(ms),plot(ms,6),tbl(ms,6)

        common /byteflip/ doflip,doflop
        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime
        common /iwin/ iwin1,iwin2

        print '(" Total effective integration time per beam:",f8.2)', &
          sumtime
        write(iupl,'(" Summed addtime:",f8.2)') sumtime
        write(iusumh,'("sumtime = ",f8.2)') sumtime
        comment = 'summed effective time (s)'
        if(dofits) call fithreal('SUMTIME ',sumtime,comment,iufish)
        write(iusumh,'("iunits  = erg/s cm2 sr cm-1")')
        write(iusumh,'("funits  = Jy")')
        write(iusumh,'("bitpix  = -32")')

      if((modeinst==MHIMED).or.(modeinst==MHILOW)) then

        ns = ny*norder
        if(ns>ms) then
          print '("*** Clipping off orders to fit in spec array",2i4)', &
            norder,ns/ny
          norder = ns/ny
          ns = norder*ny
        endif
        lt = nt+4
        if(lt>mt) then
          lt = mt
          print '("***Why was nt > mt?")'
        endif

        write(iusumh,'("nx      = ",i6)') ns
        write(iusumh,'("ny      = ",i6)') lt
        write(iusumh,'("nz      = 1")')

        line = 4*ns
        if(rdfits) then
          open(unit=iuredd,status='scratch',access='direct', &
            recl=line,err=90)
        else
          open(unit=iuredd,file=outfile,access='direct', &
            recl=line,err=90)
        endif
        do it = 1,lt
          if(doflop) then
            call flipr4(specsum(1,it),temparr,ns)
            write(unit=iuredd,rec=it,err=90) (temparr(is),is=1,ns)
          else
            write(unit=iuredd,rec=it,err=90) (specsum(is,it),is=1,ns)
          endif
        enddo   

        if(dofits) then
          bunits = 'erg/s cm2 sr cm-1'
          comment = 'cgs intensity units'
          call fithchar('BUNIT   ',bunits,comment,iufish)
          if(index(fitsfile,' ')==1) return
          call redstart(ns,nt,1,fitsfile,iufish,ierr)
          if(ierr>0) return

          do it = 5,lt
            call fitsrdat(iufred,specsum(1,it),ns)
          enddo   
          call fitsrpad(iufred)

!          do it = 1,4
!            call redxtend(iufred,ns,1,0,it)
!            call fitsrdat(iufred,specsum(1,it),ns)
!            call fitsrpad(iufred)
!          enddo   

          do is = 1,ns
            do it = 1,4
              tbl(is,it) = specsum(is,it)
            enddo
            tbl(is,5) = 1.0e+04/specsum(is,1)
          enddo
          call tblxtend(iufred,ns,5,tbl)

          call fitsclos(iufred)
        endif

        dopp = (1.+radvel/CLIGHT)
        it1 = nt/4+1
        it2 = 3*nt/4
        do is = 1,ns
          iorder = (is+ny-1)/ny
          iy = mod(is-1,ny)+1
          blksum = 0.
          isum = 0
          do it = it1,it2
            ix = (iorder-1)*nt+it
            if((ix<=nx).and.(cards(ix,iy,3)>0.)) then
              blksum = blksum+cards(ix,iy,3)
              isum = isum+1
            endif
          enddo   
          plot(is,1) = dopp*specsum(is,1)
          plot(is,2) = specsum(is,2)
!	did I put blk in plot(is,3) for atfit nlin correction?
!          plot(is,3) = blksum/max(1,isum)
          plot(is,3) = specsum(is,3)
          plot(is,4) = specsum(is,4)
          plot(is,5) = blksum/max(1,isum)
        enddo   
        scale = -1.
        if(ask) scale = 0.
        nplt = norder
        print '(" Plotting summed spectrum")'
        call plots(ms,ny,nplt,plot,scale,outfile,iwin2)
        call waitasec(pause,.false.,specsum,ms,ns,ny)
        call sleep(2)

      elseif((modeinst==MMED).or.(modeinst==MLOW)) then

        write(iusumh,'("nx      = ",i6)') nx
        write(iusumh,'("ny      = ",i6)') ny
        write(iusumh,'("nz      = 1")')

        line = 4*nx
        if(rdfits) then
          open(unit=iuredd,status='scratch',access='direct', &
            recl=line,err=90)
        else
          open(unit=iuredd,file=outfile,access='direct', &
            recl=line,err=90)
        endif
        do iy = 1,ny
          if(doflop) then
            call flipr4(arrsum(1,iy),temparr,nx)
            write(unit=iuredd,rec=iy,err=90) (temparr(ix),ix=1,nx)
          else
            write(unit=iuredd,rec=iy,err=90) (arrsum(ix,iy),ix=1,nx)
          endif
        enddo   

        if(dofits) then
          bunits = 'erg/s cm2 sr cm-1'
          comment = 'cgs intensity units'
          call fithchar('BUNIT   ',bunits,comment,iufish)
          call redstart(nx,ny-4,1,fitsfile,iufish,ierr)
          if(ierr>0) return

          do iy = 5,ny
            call fitsrdat(iufred,arrsum(1,iy),nx)
          enddo   
          call fitsrpad(iufred)

!          do iy = 1,4
!            call redxtend(iufred,nx,1,0,iy)
!            call fitsrdat(iufred,arrsum(1,iy),nx)
!            call fitsrpad(iufred)
!          enddo   

          do ix = 1,nx
            do iy = 1,4
              tbl(ix,iy) = arrsum(ix,iy)
            enddo
            tbl(ix,5) = 1.0e+04/arrsum(ix,1)
          enddo
          call tblxtend(iufred,nx,5,tbl)

          call fitsclos(iufred)
        endif

        dopp = (1.+radvel/CLIGHT)
        do ix = 1,nx
          plot(ix,1) = dopp*arrsum(ix,1)
          plot(ix,2) = arrsum(ix,2)
          plot(ix,3) = arrsum(ix,3)
          plot(ix,4) = arrsum(ix,4)
          blksum = 0.
          isum = 0
          do iy = 1,ny
            if(cards(ix,iy,3)>0.) then
              blksum = blksum+cards(ix,iy,3)
              isum = isum+1
            endif
          enddo   
          plot(ix,5) = blksum/max(1,isum)
        enddo   
        nplt = 1
        scale = -1.
        if(ask) scale = 0.
        print '(" Plotting summed spectrum")'
        call plots(ms,nx,nplt,plot,scale,outfile,iwin2)
        call waitasec(pause,.false.,arrsum,mx,nx,ny)
        call sleep(2)

      else

        write(iusumh,'("nx      = ",i6)') nx
        write(iusumh,'("ny      = ",i6)') ny
        write(iusumh,'("nz      = 1")')

        line = 4*nx
        if(rdfits) then
          open(unit=iuredd,status='scratch',access='direct', &
            recl=line,err=90)
        else
          open(unit=iuredd,file=outfile,access='direct', &
            recl=line,err=90)
        endif
        do iy = 1,ny
          if(doflop) then
            call flipr4(arrsum(1,iy),temparr,nx)
            write(unit=iuredd,rec=iy,err=90) (temparr(ix),ix=1,nx)
          else
            write(unit=iuredd,rec=iy,err=90) (arrsum(ix,iy),ix=1,nx)
          endif
        enddo   

        if(dofits) then
          bunits = 'erg/s cm2 sr cm-1'
          comment = 'cgs intensity units'
          call fithchar('BUNIT   ',bunits,comment,iufish)
          if(index(fitsfile,' ')<=1) return
          call redstart(nx,ny,1,fitsfile,iufish,ierr)
          if(ierr>0) return

          do iy = 1,ny
            call fitsrdat(iufred,arrsum(1,iy),nx)
          enddo   
          call fitsrpad(iufred)

          call fitsclos(iufred)
        endif

        if(pause.and.ask) then
          print '(" Plotting summed image")'
          call grey(mx,nx,ny,arrsum,-1.,-1.,iwin1)
          call waitasec(pause,.true.,arrsum,mx,nx,ny)
        elseif(verbose) then
          print '(" Plotting summed image")'
          call grey(mx,nx,ny,arrsum,0.,0.,iwin1)
          call waitasec(pause,.true.,arrsum,mx,nx,ny)
        endif

      endif

        ierr = 0
        close(unit=iuredd)
        close(unit=iusumh)
        return

   90   print '(" Error writing sum to ",a24)', outfile
        ierr = 9
        close(unit=iuredd)
        close(unit=iusumh)
        return

      end


      subroutine showsum(ierr)
!	plot summed spectrum

        use dims
        use modes

        real nodpa,krot,lores,kmirror
        logical baddata
        logical doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits
        logical abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause
        logical testrun,ffttort,kfix,rdhalf,irtf
        character(32) nomon

        common /nn/ nx,ny,nc,norder,ns,nt
        common /hdrf/ waveno0,temp,dist,airmass,slitpa,nodpa, &
                frtime,obstime,tottime,beamtime,addtime, &
                echelle,lores,kmirror,slit,filter,calwheel, &
                gain,pixelwd,omegap,efl,date,slitwid,pltscl
        common /flags/ modeobs,modeinst,modecard,modeblk,modetort,baddata, &
                doshift,doaddwt,dosubsky,dosum,doscansum,domon,dofits,rdfits, &
                abba,nodon,sincwt,crossdisp,rdrst,fowler,verbose,ask,pause, &
                testrun,ffttort,kfix,rdhalf,irtf
        common /torts/ slitrot,krot,detrot,spacing,xorder1,wno0,radvel, &
                hrfl,hrr,hrg,hrdgr,xdfl,xdr,xdg,xddgr,brl,x0brl,y0brl
        common /ss/ arrsum(mx,my),specsum(ms,mt),scansum(mx,my,mp), &
                cards(mx,my,4),sumtime
        common /iwin/ iwin1,iwin2

      nomon = ' '

      if((modeinst==MHIMED).or.(modeinst==MHILOW)) then

        scale = -1.
        if(ask) scale = 0.
        nplt = norder
        print '(" Plotting summed spectrum")'
        call plots(ms,ny,nplt,specsum,scale,nomon,iwin2)
        call waitasec(pause,.false.,specsum,ms,ns,ny)

      elseif((modeinst==MMED).or.(modeinst==MLOW)) then

        nplt = 1
        scale = -1.
        if(ask) scale = 0.
        print '(" Plotting summed spectrum")'
        call plots(ms,nx,nplt,arrsum,scale,nomon,iwin2)
        call waitasec(pause,.false.,arrsum,mx,nx,ny)

      else

        print '(" Plotting summed image")'
        if(pause.and.ask) then
          call grey(mx,nx,ny,arrsum,-1.,-1.,iwin1)
        else
          call grey(mx,nx,ny,arrsum,0.,0.,iwin1)
        endif
        call waitasec(pause,.true.,arrsum,mx,nx,ny)

      endif

        ierr = 0
        return

      end


      function bnu(w,t,ierr)
!	black-body function

        use consts

        bnu = 2.*hc2*w**3/(sexp(hck*w/t)-1.)

        return
      end


      function pnu(w,t,ierr)
!	black-body photon function

        use consts

        pnu = 2.*(CLIGHT*1.0e+05)*w**2/(sexp(hck*w/t)-1.)

        return
      end


      subroutine smooth(nx,ny,arr,sarr,sx,sy)
!	apparently not used

        use dims

        dimension arr(mx,my),sarr(mx,my)
        dimension aline(mx),sline(mx)

        if(sx<=1.) goto 201
        xnx = nx-1
        do iy = 1,ny
          do ix = 1,nx
            aline(ix) = arr(ix,iy)
          enddo   
          sline(1) = aline(1)
          sline(nx) = aline(nx)
          do ix = 2,nx-1
            x = ix-1
            sx2 = min((sx-1.)/2.,x,xnx-x)
            nsx = sx2
            sum = (sx2-nsx)*(aline(ix-nsx-1)+aline(ix+nsx+1))
            do isx = ix-nsx,ix+nsx
              sum = sum+aline(isx)
            enddo   
            sline(ix) = sum/sx
          enddo   
          do ix = 1,nx
            sarr(ix,iy) = sline(ix)
          enddo   
        enddo   

  201   if(sy<=1.) goto 301
        yny = ny-1
        do ix = 1,nx
          do iy = 1,ny
            aline(iy) = sarr(ix,iy)
          enddo   
          sline(1) = aline(1)
          sline(ny) = aline(ny)
          do iy = 2,ny-1
            y = iy-1
            sy2 = min((sy-1.)/2.,y,yny-y)
            nsy = sy2
            sum = (sy2-nsy)*(aline(iy-nsy-1)+aline(iy+nsy+1))
            do isy = iy-nsy,iy+nsy
              sum = sum+aline(isy)
            enddo   
            sline(iy) = sum/sy
          enddo   
          do iy = 1,ny
            sarr(ix,iy) = sline(iy)
          enddo   
        enddo   

  301   return

      end


      function ibad(val)
!	check for NaN or Inf

        if(val/=val) then
          ibad = 1
        elseif((val/=0.).and.(val/val/=1.)) then
          ibad = 2
        else
          ibad = 0
        endif

        return
      end


      subroutine checkarr(arr,name,ierr)
!	check array for NaN or Inf

        use dims

        real arr(mx,my)
        real plot(mx,my)
        character(8) name
        character(8) yn
        logical good,setillum

        common /nn/ nx,ny,nc,norder,ns,nt
        common /goods/ good(mx,my),illum(mx,my),illx(mx),wt(mp)
        common /iwin/ iwin1,iwin2

        setillum = (ierr==0)
        ierr = 0
        jerr = 0
        do iy = 1,ny
          do ix = 1,nx
            plot(ix,iy) = illum(ix,iy)
            if(ibad(arr(ix,iy))>0) then
              if(ierr<4) print '("*** ",a8,"(",2i4,") = ",e10.2)', &
                name,ix,iy,arr(ix,iy)
              arr(ix,iy) = 0.
              illum(ix,iy) = 0
              plot(ix,iy) = -2.
              ierr = ierr+1
            elseif((arr(ix,iy)==0.).and.(illum(ix,iy)==1)) then
              if(setillum) illum(ix,iy) = 0
              plot(ix,iy) = 2.
              jerr = jerr+1
            endif
          enddo   
        enddo   
        if(ierr>0) then
          print '(i6," bad pixels in ",a8)', ierr,name
          call grey(mx,nx,my,plot,0.,0.,iwin1)
          call waitasec(.true.,.true.,plot,mx,nx,ny)
        endif
        if((.not.setillum).and.(jerr>0)) then
          print '(i6," pixels = 0 in ",a8, &
            " where illum = 1")', jerr,name
        endif

        ierr = 0
        return
      end


      subroutine matinv(arr,ni,det)
!	invert a symmetric array

        integer, parameter :: mi = 256

        real arr(ni,ni)
        integer ik(mi),jk(mi)

        if(ni>mi) then
          print '(" Error in matinv.  Array too big.")'
          det = 0.
          return
        endif

        det = 1.
        do k = 1,ni
          armax = 0.
          abmax = 0.
  102     do i = k,ni
            do j = k,ni
              if(abs(arr(i,j))>abmax) then
                armax = arr(i,j)
                abmax = abs(armax)
                ik(k) = i
                jk(k) = j
              endif
            enddo   
          enddo   
          if(abmax==0.) then
            print '(" Error in matinv.  Matrix singular.")'
            det = 0.
            return
          endif

          i = ik(k)
          if(i<k) then
            goto 102
          elseif(i>k) then
            do j = 1,ni
              save = arr(k,j)
              arr(k,j) = arr(i,j)
              arr(i,j) = -save
            enddo   
          endif
          j = jk(k)
          if(j<k) then
            goto 102
          elseif(j>k) then
            do i = 1,ni
              save = arr(i,k)
              arr(i,k) = arr(i,j)
              arr(i,j) = -save
            enddo   
          endif

          do i = 1,ni
            if(i/=k) arr(i,k) = -arr(i,k)/armax
          enddo   
          do i = 1,ni
            do j = 1,ni
              if((i/=k).and.(j/=k)) &
                arr(i,j) = arr(i,j)+arr(i,k)*arr(k,j)
            enddo   
          enddo   
          do j = 1,ni
            if(j/=k) arr(k,j) = arr(k,j)/armax
          enddo   
          arr(k,k) = 1./armax
          det = det*armax

        enddo   

        do l = 1,ni
          k = ni-l+1
          j = ik(k)
          if(j>k) then
            do i = 1,ni
              save = arr(i,k)
              arr(i,k) = -arr(i,j)
              arr(i,j) = save
            enddo   
          endif
          i = jk(k)
          if(i>k) then
            do j = 1,ni
              save = arr(k,j)
              arr(k,j) = -arr(i,j)
              arr(i,j) = save
            enddo   
          endif
        enddo   

        do i = 1,ni
          do j = 1,ni
            if(arr(i,j)/=arr(i,j)) then
              print '(" Error in matinv.  Probable overflow.")'
              det = 0.
              return
            endif
          enddo   
        enddo   

        return
      end

