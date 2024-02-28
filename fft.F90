      subroutine fft(fr,fi,nu)

        common /trigtab/costab(0:4095),sintab(0:4095)
        real fr(2**nu),fi(2**nu)
        save ntt

        n = 2**nu
        if(n/=ntt) then
          call calctt(n)
          ntt = n
        endif

        k = 0
        n2 = n/2
        do l = 1,nu
 100      do i = 1,n2
            itab = ibitr(k/n2,nu)
            c = costab(itab)
            s = sintab(itab)
            k = k+1
            k2 = k+n2
            tr = fr(k2)*c+fi(k2)*s
            ti = fi(k2)*c-fr(k2)*s
            fr(k2) = fr(k)-tr
            fi(k2) = fi(k)-ti
            fr(k) = fr(k)+tr
            fi(k) = fi(k)+ti
          enddo
          k = k+n2
          if(k<n) go to 100
          k = 0
          n2 = n2/2
        enddo

        do k = 1,n
          i = ibitr(k-1,nu)+1
          if(i<=k) cycle
          tr = fr(k)
          ti = fi(k)
          fr(k) = fr(i)
          fi(k) = fi(i)
          fr(i) = tr
          fi(i) = ti
        enddo

        return
      end


        function ibitr(j,nu)

        j1 = j
        ibitr = 0

        do i = 1,nu
          j2 = j1/2
          ibitr = ibitr*2+(j1-2*j2)
          j1 = j2
        enddo

        return
      end


      subroutine calctt(n)

        common /trigtab/costab(0:4095),sintab(0:4095)

        twopin = 6.283185/n
        costab(0) = 1.
        sintab(0) = 0.
        do i = 1,n/2
          arg = twopin*i
          costab(i) = cos(arg)
          costab(n-i) = costab(i)
          sintab(i) = sin(arg)
          sintab(n-i) = -sintab(i)
        enddo

        return
      end


      subroutine ifft(fr,fi,nu)

        real fr(2**nu),fi(2**nu)

        n = 2**nu
        do i = 1,n
          fi(i) = -fi(i)
        enddo
        call fft(fr,fi,nu)
        do i = 1,n
          fr(i) = fr(i)/n
          fi(i) = -fi(i)/n
        enddo

        return
      end


      subroutine rfft(fr,fi,nu)

        real fr(2**nu),fi(2**nu)

        n = 2**nu
        mu = nu-1
        m = 2**mu
        pim = 3.141593/m

        do j = 1,m
          i = 2*j-1
          fr(j) = fr(i)
          fi(j) = fr(i+1)
        enddo

        call fft(fr,fi,mu)

        do j = 2,m
          j1 = m-j+2
          i1 = n-j+2
          c = cos(pim*(j-1))
          s = sin(pim*(j-1))
          fr(i1) = ((fr(j)+fr(j1))+c*(fi(j)+fi(j1))-s*(fr(j)-fr(j1)))/2.
          fi(i1) = -((fi(j)-fi(j1))-s*(fi(j)+fi(j1))-c*(fr(j)-fr(j1)))/2.
        enddo

        do i = 2,m
          i1 = n-i+2
          fr(i) = fr(i1)
          fi(i) = -fi(i1)
        enddo
        fr(m+1) = fr(1)-fi(1)
        fi(m+1) = 0.
        fr(1) = fr(1)+fi(1)
        fi(1) = 0.

        return
      end


      subroutine tdfft(gr,gi,nux,nuy)

        real gr(2**(nux+nuy)),gi(2**(nux+nuy))
        real fr(4096),fi(4096)

        nx = 2**nux
        ny = 2**nuy

        do ix = 1,nx
          do iy = 1,ny
            fr(iy) = gr(nx*(iy-1)+ix)
            fi(iy) = gi(nx*(iy-1)+ix)
          enddo

          call fft(fr,fi,nuy)

          do iy = 1,ny
            gr(nx*(iy-1)+ix) = fr(iy)
            gi(nx*(iy-1)+ix) = fi(iy)
          enddo
        enddo

        do iy = 1,ny
          i1 = nx*(iy-1)+1
          call fft(gr(i1),gi(i1),nux)
        enddo

        return
      end


      subroutine tdifft(gr,gi,nux,nuy)

        real gr(2**(nux+nuy)),gi(2**(nux+nuy))
        real fr(4096),fi(4096)

        nx = 2**nux
        ny = 2**nuy
        xn = nx*ny

        do ix = 1,nx
          do iy = 1,ny
            fr(iy) = gr(nx*(iy-1)+ix)
            fi(iy) = -gi(nx*(iy-1)+ix)
          enddo

          call fft(fr,fi,nuy)

          do iy = 1,ny
            gr(nx*(iy-1)+ix) = fr(iy)
            gi(nx*(iy-1)+ix) = fi(iy)
          enddo
        enddo

        do iy = 1,ny
          i1 = nx*(iy-1)+1

          call fft(gr(i1),gi(i1),nux)

          do ix = 1,nx
            ii = nx*(iy-1)+ix
            gr(ii) = gr(ii)/xn
            gi(ii) = -gi(ii)/xn
          enddo
        enddo

        return
      end


      subroutine tdrfft(gr,gi,nux,nuy)

        real gr(2**(nux+nuy)),gi(2**(nux+nuy))
        real fr(4096),fi(4096)

        nx = 2**nux
        ny = 2**nuy
        ny2 = ny/2

!	real fft in y for each x

        do ix = 1,nx
          do iy = 1,ny
            fr(iy) = gr(nx*(iy-1)+ix)
            fi(iy) = 0.
          enddo

          call rfft(fr,fi,nuy)

          do iy = 1,ny
            gr(nx*(iy-1)+ix) = fr(iy)
            gi(nx*(iy-1)+ix) = fi(iy)
          enddo
        enddo

!	fft in x for each y

!	g(iy=1), g(iy=ny2+1) are real

        call rfft(gr(1),gi(1),nux)
        call rfft(gr(ny2*nx+1),gi(ny2*nx+1),nux)

!	g(-iy) = g*(iy), so G(-iy) = G*(iy,reversed)

        do iy=2,ny2
          i1 = nx*(iy-1)+1
          call fft(gr(i1),gi(i1),nux)
!	zeroth point of row iy
          i1 = nx*(iy-1)
!	last point of row -iy
          i2 = nx*(ny-iy+2)
          gr(i2-nx+1) = gr(i1+1)
          gi(i2-nx+1) = -gi(i1+1)
          do ix = 2,nx
            gr(i2-ix+2) = gr(i1+ix)
            gi(i2-ix+2) = -gi(i1+ix)
          enddo
        enddo

        return
      end


      subroutine rpow(fr,nu)

        real fr(2**nu)
        real fi(2**nu)

        n = 2**nu
        n2 = n/2

        call rfft(fr,fi,nu)

        do i = 2,n2-1
          fr(i) = fr(i)**2+fi(i)**2+fr(n-i)**2+fi(n-i)**2
          fr(n2+i) = 0.
        enddo
        fr(1) = fr(1)**2+fi(1)**2
        fr(n2) = 0.
        fr(n2+1) = 0.
        fr(n) = 0.

        return
      end


      function sexp(arg)
!	safe and relatively efficient exp function
 
        parameter (big = 72.)
        parameter (small = 0.01)
 
        if(arg<(-big)) then
          sexp = 1.0e-31
        elseif(arg>big) then
          sexp = 1.0e+31
        elseif(abs(arg)<small) then
          sexp = 1.0+arg+arg**2/2.+arg**3/6.
        elseif((arg<0.).or.(arg>0.)) then
          sexp = exp(arg)
        else
          print '(es10.2," passed to sexp")', arg
          sexp = 0.
        endif
 
        return
      end
 
 
      function slog(arg)
!	safe and efficient log function
 
        parameter (small = 0.01)
        parameter (huge = 1.0e+31)
        parameter (tiny = 1.0e-32)
 
        if(arg<=0.) then
          slog = 0.
        elseif(abs(arg-1.)<small) then
          slog = arg-1.
          slog = slog-slog**2/2.+slog**3/3.
        elseif(arg>0.) then
          slog = alog(arg)
        else
          print '(es10.2," passed to slog")', arg
          slog = 0.
        endif
 
        return
      end


