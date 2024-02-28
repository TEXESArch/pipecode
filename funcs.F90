      function bw(w,t)
!	black body function

        parameter (hck = 1.439)
        parameter (hc2 = 5.954e-06)

        arg = hck*w/t
        if(arg<=0.) then
          bw = 0.
        elseif(arg<100.) then
          bw = 2.*hc2*w**3/(texp(arg)-1.)
        else
          bw = 0.
        endif

      end function bw


      function texp(arg)
!	exp(arg) using look-up table

        parameter (hun = 100.)
        common /exptab/ expi(10239),too

        if(abs(arg)<hun) then
          harg = abs(hun*arg)
          iarg = int(harg)
          darg = (harg-iarg)/hun
          texp = (1.-darg+darg**2/too)*expi(iarg+1)
          if(arg>0.) texp = 1./texp
        elseif(arg<0.) then
          texp = 0.
        else
          texp = 1./expi(10000)
        endif

        return
 
      entry expinit

        do iarg = 1,10000
          expi(iarg) = sexp(-(iarg-1)/hun)
        enddo
        darg = 1./hun
        too = darg**2/(exp(-darg)-1.+darg)
        expinit = 1.

      end function texp


!      function sexp(arg)
!!	safe and relatively efficient exp function
!
!        parameter (hun = 70.)
!        parameter (small = 0.01)
!        parameter (huge = 0.999e+31)
!        save nan
!        if((nan<0).or.(nan>8)) nan = 0
!
!        if(arg<(-hun)) then
!          sexp = 0.
!        elseif(arg>hun) then
!          sexp = huge
!        elseif(abs(arg)<small) then
!          sexp = 1.+arg+arg**2/2.+arg**3/6.
!        elseif((arg<0.).or.(arg>0.)) then
!          sexp = exp(arg)
!        else  !  NaN
!          print '(e10.2," passed to sexp",i6)', arg,nan
!          nan = nan+1
!          if(nan>8) stop
!          sexp = -1.
!        endif
!
!      end function sexp


!      function slog(arg)
!!	safe log function
!
!        parameter (huge = 0.999e+31)
!        parameter (tiny = 1.0e-32)
!        parameter (small = 0.01)
!
!        save nan
!        if((nan<0).or.(nan>8)) nan = 0
!
!        if(arg<=0.) then
!          slog = -100.
!        elseif(arg<tiny) then
!          slog = -72.
!        elseif(arg>huge) then
!          slog = 72.
!        elseif(abs(arg-1.)<small) then
!          slog = arg-1.
!          slog = slog-slog**2/2.+slog**3/3.
!        elseif(arg>tiny) then
!          slog = log(arg)
!        else  !  NaN
!          print '(e10.2," passed to slog",i6)', arg
!          nan = nan+1
!          if(nan>8) stop
!          slog = 0.
!        endif
!
!      end function slog


      function slog10(arg)
!	safe log_10 function

        parameter (huge = 0.999e+31)
        parameter (tiny = 1.0e-32)
        parameter (small = 0.01)
        parameter (aloge = 0.434294)

        save nan
        if((nan<0).or.(nan>8)) nan = 0

        if(arg<=0.) then
          slog10 = -32.
        elseif(arg<tiny) then
          slog10 = -32.
        elseif(arg>huge) then
          slog10 = 31.
        elseif(abs(arg-1.)<small) then
          slog = arg-1.
          slog10 = aloge*(slog-slog**2/2.+slog**3/3.)
        elseif(arg>tiny) then
          slog10 = log10(arg)
        else  !  NaN
          print '(e10.2," passed to slog10",i6)', arg
          nan = nan+1
          if(nan>8) stop
          slog10 = 0.
        endif

      end function slog10


      function satan(arg)
!	safe atan function

        parameter (small = 0.01)
        parameter (big = 100.)
        parameter (pi2 = 1.570796)

        save nan
        if((nan<0).or.(nan>8)) nan = 0

        if(abs(arg)<small) then
          satan = arg-arg**3/3.+arg**5/5.
        elseif(abs(arg)>big) then
          arginv = abs(1./arg)
          satan = pi2-arginv+arginv**3/3.-arginv**5/5.
          if(arg<0.) satan = -satan
        elseif(abs(arg)>0.) then
          satan = atan(arg)
        else  !  NaN
          print '(e10.2," passed to satan",i6)', arg
          nan = nan+1
          if(nan>8) stop
          satan = 0.
        endif

      end function satan


      logical function badnum(val)
!	test val for NaN or Inf

        if(val==0.) then
          badnum = .false.
        else
          test = val/val
          badnum = (test/=test)
        endif

      end function badnum


      function quartf(arr,nx,iq,ierr)
!	return iqth quartile of arr

        real arr(nx),arrf(65536)
        logical badnum

        nxq = nx
        arrf(1) = arr(1)
        if(arrf(1)==0.) nxq = nxq-1
        do ix = 1,nx
          arri = arr(ix)
          if(badnum(arri)) then
            ierr = 1
            arri = 0.
          endif
          if(arri==0.) then
            nxq = nxq-1
            arrf(ix) = arri
            goto 100
          endif
          do jx = ix,2,-1
            arrj = arrf(jx-1)
            if((arrj==0.).or.(arri<arrj)) then
              arrf(jx) = arrj
            else
              arrf(jx) = arri
              goto 100
            endif
          enddo
          arrf(1) = arri
  100   enddo
        if(iq==0) then
          quartf = arrf(1)
        elseif(iq==4) then
          ixq = nxq-(nxq/64)
          quartf = arrf(ixq)
        else
          ixq = (iq*nxq+3)/4
          if((mod(nxq,4)==0).or.(iq==2).or.(mod(nxq,2)==0)) then
            quartf = (arrf(ixq)+arrf(ixq+1))/2.
          else
            quartf = arrf(ixq)
          endif
        endif

      end function quartf


      subroutine quartz(arr2,arr3,nx,ny,nz,iq,ierr)
!	collapse arr3 in z by quartile filtering

        use dims

        real arr2(nx,ny),arr3(nx,ny,nz)
        real arrz(mp)

        ierr = 0
        do iy = 1,ny
          do ix = 1,nx
            do iz = 1,nz
              arrz(iz) = arr3(ix,iy,iz)
            enddo
            arr2(ix,iy) = quartf(arrz,nz,iq,ierr)
          enddo
        enddo

      end subroutine quartz


