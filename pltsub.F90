!
!	Note: to enable screen plotting, calling program must have:
!
!	integer pgopen
!	idev = pgopen('[file]/type')
!	if(idev<=0) complain
!	...
!	call pgclos
!


      subroutine grey(nx,lx,ny,z,zmn,zmx,iwin)
!	zmn = zmx: linear min->max
!	zmn < zmx: linear zmn->zmx
!	zmn=0 zmx<0: sqrt 0->|zmx|
!	zmn=0 zmx=-1: sqrt min->max
!	zmn > zmx: log min->max

        integer, parameter :: mx=2048
        integer, parameter :: my=2048
        dimension z(1)
        dimension a(mx,my),tr(6)
        character*8 notval
        logical dosqrt,dolog

        save xannot,yannot

        if(iwin<1) return
        call pgslct(iwin)
        call pgsch(0.)
        call pgsubp(1,1)

        ll = max0(lx,ny)
        if((nx>mx).or.(ll>my)) then
          print '(" Image too big for grey",5i6)', nx,lx,ny,mx,my
          return
        endif
        tr(1) = -0.5/ll
        tr(2) = 1./ll
        tr(3) = 0.
        tr(4) = -0.5/ll
        tr(5) = 0.
        tr(6) = 1./ll
        xannot = (0.9*lx)/ll
        yannot = (0.9*ny)/ll

        zmax = -1.0e+06
        zmin = 1.0e+06
        iz = 0
        do iy=1,my
          do ix=1,mx
            a(ix,iy) = 0.
            if((ix<=nx).and.(iy<=ny)) then
              iz = iz+1
              zi = z(iz)
              if(zi>zmax) zmax = zi
              if(zi<zmin) zmin = zi
            endif
          enddo
        enddo

        dosqrt = ((zmn==0.).and.(zmx<0.))
        dolog = ((zmn.ne.0.).and.(zmx<zmn))
        iz = 0
        do iy=1,ny
          do ix=1,nx
            iz = iz+1
            zi = z(iz)
            if(dosqrt) then
              a(ix,iy) = sqrt(zi-zmin)
            elseif(dolog) then
              a(ix,iy) = sign(-zmx*alog(1.+abs(zi/zmx)),zi)
            else
              a(ix,iy) = zi
            endif
          enddo
        enddo

        if(dosqrt) then
          if(zmx==-1.) then
            zmax = 0.
            zmin = 0.
          else
            zmax = sqrt(abs(zmx-zmin))
            zmin = 0.
          endif
        elseif(dolog) then
          zmin = 0.
          zmax = 0.
        else
          zmin = zmn
          zmax = zmx
        endif

   24   if(zmin==zmax) then
          do iy=2,ny-1
            do ix=2,lx-1
              if(a(ix,iy).ne.a(ix,iy)) cycle
              if((a(ix-1,iy)<zmin).and.(a(ix,iy)<zmin) &
                .and.(a(ix+1,iy)<zmin)) &
                zmin = max(a(ix-1,iy),a(ix,iy),a(ix+1,iy))
              if((a(ix-1,iy)>zmax).and.(a(ix,iy)>zmax) &
                .and.(a(ix+1,iy)>zmax)) &
                zmax = min(a(ix-1,iy),a(ix,iy),a(ix+1,iy))
            enddo
          enddo
        elseif(zmax==0.) then
          do iy=2,ny-1
            do ix=2,lx-1
              if((a(ix-1,iy)>zmax).and.(a(ix,iy)>zmax) &
                .and.(a(ix+1,iy)>zmax)) &
                zmax = min(a(ix-1,iy),a(ix,iy),a(ix+1,iy))
            enddo
          enddo
        elseif(zmax<zmin) then
          print '(" Enter zmin, zmax for grayscale plot")'
          print '(" 0,0 to get both from image; x,0 to get only zmax")'
          read *, zmin, zmax
          dosqrt = (zmax<0.)
          if(zmax.ne.0.) then
            zmn = zmin
            zmx = zmax
            if(dosqrt) zmax = abs(zmax)
          endif
          goto 24
        endif
        if((0.*zmin.ne.0.).or.(0.*zmax.ne.0.)) &
          print '(" In grey zmin,zmax =",1p,2e10.2)', zmin,zmax

        if(zmax>zmin) then
          call pgenv(0.,1.,0.,1.,1,-2)
          call pggray(a,mx,my,1,ll,1,ll,zmax,zmin,tr)
        else
          print '("*** In grey: zmin, zmax = ",1p,2e10.2)', zmin,zmax
        endif

        return

      entry annot(val)

          write (unit=notval,fmt='(f8.2)') val
!          call pgtext(xannot,yannot,notval)
!          call pgptxt(xannot,yannot,-90.,0.,notval)
          call pgptxt(xannot,yannot,0.,0.,notval)

        return

      end



      subroutine contour(nx,lx,ny,z,dz0,iwin)

        character*8 yn
        character*24 monfile
        character*32 levels
        logical ask,start,transp,mirx,miry,sqt,peak
        dimension z(1)
        dimension zsum(1024)

        real, parameter :: alog2=0.30103

        if(iwin<1) return
        call pgslct(iwin)
        call pgsch(1.)
        call pgsubp(1,1)

        ask = (dz0==0.)

        npt = nx*ny
        ixmn = 1
        ixmx = lx-1
        iymn = 1
        iymx = ny-1

  10    zmin = z(1)
        zmax = z(1)
        do ipt=1,npt
          if(z(ipt).ne.z(ipt)) then
            z(ipt) = 0.
            cycle
          endif
          zmin = min(z(ipt),zmin)
          zmax = max(z(ipt),zmax)
        enddo
        if(zmax==zmin) then
          print '(" zmax = zmin = ",1p,e10.2)', zmax
          return
        endif
        if(dz0<=0.) then
          al10 = log10(zmax-zmin)
          il10 = int(al10)
          al10 = al10-il10
          il2 = int(al10/alog2)
          dz = (2.**il2)*(10.**(il10-1))
        else
          dz = dz0
        endif
        nz = (zmax-zmin)/dz+1.
        if(zmin<0.) zmin = zmin-dz
        zmin = dz*aint(zmin/dz)

        if(ask) then
          print '(/" xmn, xmx, ymn, ymx =",4i8," OK? (y)")', &
            ixmn,ixmx+1,iymn,iymx+1
          read '(a8)', yn
          if(index(yn,'n')>0) then
  12        print '(" Enter xmn, xmx, ymn, ymx")'
            read *, ixmn,ixmx,iymn,iymx
            if((ixmn<1).or.(ixmx>nx) &
              .or.(iymn<1).or.(iymx>ny)) then
              print '(" Values out of range.")'
              go to 12
            endif
            ixmx = ixmx-1
            iymx = iymx-1
          endif
          print '(" Mirror or transpose plot? (n)")'
          read '(a8)', yn
          if(index(yn,'y')>0) then
            print '(" Enter trans, xmir, ymir (T = 1)")'
            read *, trans,xmir,ymir
            mirx = (xmir==1.)
            miry = (ymir==1.)
            xymir = 1.
            if(mirx) xymir = -xymir
            if(miry) xymir = -xymir
            transp = (trans==1.)
          else
            mirx = .false.
            miry = .false.
            xymir = 1.
            transp = .false.
          endif
        endif
        if(transp) then
          xmn = iymn-0.5
          xmx = iymx+1.5
          ymn = ixmn-0.5
          ymx = ixmx+1.5
        else
          xmn = ixmn-0.5
          xmx = ixmx+1.5
          ymn = iymn-0.5
          ymx = iymx+1.5
        endif
        if(mirx) then
          smn = xmn
          xmn = -xmx
          xmx = -smn
        endif
        if(miry) then
          smn = ymn
          ymn = -ymx
          ymx = -smn
        endif

        if(ask) then
          print '(" zmin, dz, nz =",2f10.3,i8" OK? (y)")', zmin,dz,nz
          read '(a32)', levels
  14      if(index(levels,'n')>0) then
            print '(" Enter zmin, dz, nz", &
              " (0 to compute dz or nz, dz<0 to plot sqrt)")'
            read *, zmin,dz,nz
          elseif((index(levels,'y')==0).and.(index(levels,' ')>1)) &
            then
            read(unit=levels,fmt=*) zmin,dz,nz
          endif
          sqt = (dz<0.)
          if(sqt) then
            dz = sqrt(abs(dz))
            zmin = sqrt(max(0.,zmin))
            zmax = sqrt(max(0.,zmax))
          endif
          if(dz==0.) then
            if(nz==0) go to 14
            dz = (zmax-zmin)/(nz-1)
            print '(" zmin, dz, nz = ",2f10.3,i8)', zmin,dz,nz
          elseif(nz==0) then
            nz = (zmax-zmin)/dz+1.
            print '(" zmin, dz, nz = ",2f10.3,i8)', zmin,dz,nz
          endif
        endif
        nzabs = abs(nz)
        nz2 = nzabs/2

        zmn = 0.
        zmx = 0.
        if(transp) then
          do iy = iymn,iymx
            zsumi = 0.
            do ix = ixmn,ixmx
              zsumi = zsumi+z(nx*(iy-1)+ix)
            enddo
            zsum(iy) = zsumi
            zmn = min(zmn,zsumi)
            zmx = max(zmx,zsumi)
          enddo
        else
          do ix = ixmn,ixmx
            zsumi = 0.
            do iy = iymn,iymx
              zsumi = zsumi+z(nx*(iy-1)+ix)
            enddo
            zsum(ix) = zsumi
            zmn = min(zmn,zsumi)
            zmx = max(zmx,zsumi)
          enddo
        endif

        call pgenv(xmn,xmx,(ymn-(ymx-ymn)/4.),ymx,1,0)
        call pgsls(1)

        do iy = iymn,iymx
          do ix = ixmn,ixmx
            ipt1 = nx*(iy-1)+ix
            ipt2 = nx*iy+ix
            ipt3 = nx*iy+ix+1
            ipt4 = nx*(iy-1)+ix+1
            z1 = z(ipt1)
            z2 = z(ipt2)
            z3 = z(ipt3)
            z4 = z(ipt4)
            if(sqt) then
              z1 = sqrt(max(0.,z1))
              z2 = sqrt(max(0.,z2))
              z3 = sqrt(max(0.,z3))
              z4 = sqrt(max(0.,z4))
            endif
            if((z1*z2*z3*z4)==0.) cycle
            zhi = max(max(max(z1,z2),z3),z4)
            zlo = min(min(min(z1,z2),z3),z4)
            izlo = int((zlo-zmin)/dz+2.)-1
            izhi = int((zhi-zmin)/dz+1.)-1
            izlo = max(izlo,1)
            if(nz>0) then
              izhi = min(izhi,nz)
            else
              mzlo = 0
              do
                if(izlo>=nzabs) exit
                izlo = izlo/2
                mzlo = mzlo+1
              enddo
              izlo = izlo+mzlo*nzabs/2
              mzhi = 0
              do
                if(izhi>=nzabs) exit
                izhi = izhi/2
                mzhi = mzhi+1
              enddo
              izhi = izhi+mzhi*nzabs/2
            endif
            ncon = izhi-izlo+1
            if(ncon<1) cycle
            zc = (z1+z2+z3+z4)/4.
            if(transp) then
              xa = iy
              ya = ix
              dx = 0.
              dy = 1.
            else
              xa = ix
              ya = iy
              dx = 1.
              dy = 0.
            endif
            xc = xa+0.5
            yc = ya+0.5
            if(mirx) then
              xa = -xa
              xc = -xc
              dx = -dx
            endif
            if(miry) then
              ya = -ya
              yc = -yc
              dy = -dy
            endif
            ipta = ipt1
            za = z1
            jpt = 1
            do ia = 1,2
              do ib = 1,2
                iptb = ipta+jpt
                xb = xa+dx
                yb = ya+dy
                zb = z(iptb)
                if(sqt) zb = sqrt(max(0.,zb))
                zlo = min(min(za,zb),zc)
                zhi = max(max(za,zb),zc)
                izlo = int((zlo-zmin)/dz+2.)-1
                izhi = int((zhi-zmin)/dz+1.)-1
                if(nz>0) then
                  izlo = max(izlo,1)
                  izhi = min(izhi,nz)
                else
                  mzlo = 0
                  do
                    if(izlo>=nzabs) exit
                    izlo = izlo/2
                    mzlo = mzlo+1
                  enddo
                  izlo = izlo+mzlo*nzabs/2
                  mzhi = 0
                  do
                    if(izhi>=nzabs) exit
                    izhi = izhi/2
                    mzhi = mzhi+1
                  enddo
                  izhi = izhi+mzhi*nzabs/2
                  izlo = max(izlo,1)
                endif
                do iz = izlo,izhi
                  start = .false.
                  if((nz>0).or.(iz<=nzabs)) then
                    zi = zmin+iz*dz
                  else
                    mz = (2*iz)/nzabs
                    zi = zmin+(iz-(mz-1)*nzabs/2)*(2**(mz-1))*dz
                  endif
                  if(za.ne.zb) then
                    dab = (zi-za)/(zb-za)
                    if ((dab>=0.).and.(dab<1.)) then
                      start = .true.
                      x = xa+dab*(xb-xa)
                      y = ya+dab*(yb-ya)
                      call pgmove(x,y)
                    endif
                  endif
                  if(za.ne.zc) then
                    dac = (zi-za)/(zc-za)
                    if((dac>=0.).and.(dac<1.)) then
                      x = xa+dac*(xc-xa)
                      y = ya+dac*(yc-ya)
                      if(start) then
                        call pgdraw(x,y)
                        cycle
                      else
                        call pgmove(x,y)
                      endif
                    endif
                  endif
                  if(zb.ne.zc) then
                    dbc = (zi-zb)/(zc-zb)
                    if((dbc>=0.).and.(dbc<1.)) then
                      x = xb+dbc*(xc-xb)
                      y = yb+dbc*(yc-yb)
                      call pgdraw(x,y)
                    endif
                  endif
                enddo
                jpt = nx*jpt
                sdy = dy
                dy = xymir*dx
                dx = xymir*sdy
              enddo
              ipta = ipt3
              za = z(ipta)
              if(sqt) za = sqrt(max(0.,za))
              if(mirx) then
                xa = xa-1.
              else
                xa = xa+1.
              endif
              if(miry) then
                ya = ya-1.
              else
                ya = ya+1.
              endif
              dx = -dx
              dy = -dy
              jpt = -1
            enddo
          enddo
        enddo
        if(zmn>=0.) then
          y0 = ymn-(ymx-ymn)/4.
          dydz = (ymx-ymn)/(4.*zmx)
        else
          y0 = ymn-(ymx-ymn)/8.
          dydz = (ymx-ymn)/(4.*(zmx-zmn))
        endif
        if(transp) then
          x = iymn
          y = y0+dydz*zsum(iymn)
          call pgmove(x,y)
          do ix = iymn,iymx
            x = ix
            if(mirx) x = -x
            y = y0+dydz*zsum(ix)
            call pgdraw(x,y)
          enddo
        else
          x = iabs(ixmn)
          y = y0+dydz*zsum(ixmn)
          call pgmove(x,y)
          do ix = ixmn,ixmx
            x = ix
            if(mirx) x = -x
            y = y0+dydz*zsum(ix)
            call pgdraw(x,y)
          enddo
        endif

        if(.not.ask) return

        print '(" Replot map? (n)")'
        read '(a8)', yn
        if(index(yn,'y')>0) go to 10

   39   print '(" Do you want a mongo plot file? (n)")'
        read '(a8)', yn
        if(index(yn,'y')==0) return
        print '(" Enter name of mongo file")'
        read '(a16)', monfile
        open(unit=4,file=monfile,status='new',access='sequential',err=39)
        write(4,'(a5)') 'reset'
        write(4,'(a5)') 'erase'
        write(4,'(a10)') 'expand 1.0'
        write(4,'(a9)') 'psfmode 2'
        write(4,'(a7,a16)') 'tlabel ',monfile
        write(4,'(a12)') 'window 1 1 1'
        write(4,'(a8,1p,2e12.3)') 'xlimits ',xmn,xmx
        write(4,'(a8,1p,2e12.3)') 'ylimits ',ymn,ymx
        write(4,'(a10)') 'expand 1.0'
        write(4,'(a14)') 'square 1 2 0 0'
        write(4,'(a32,3f8.3)') 'xlabel minval, interval, maxval:', &
          zmin,dz,zmax
        write(4,'(a9)') 'ptype 4 1'
        write(4,'(a8)') 'angle 45'
  
        do iy = iymn,iymx
          do ix = ixmn,ixmx
            jymn = max(1,iy-3)
            jymx = min(ny,iy+3)
            jxmn = max(1,ix-3)
            jxmx = min(lx,ix+3)
            zi = z(nx*iy+ix)
            peak = (zi>zmin+2.*abs(dz))
            if(.not.peak) jymx = jymn-1
            do jy = jymn,jymx
              do jx = jxmn,jxmx
                zj = z(nx*jy+jx)
                if(zj>zi) peak = .false.
              enddo
            enddo
            if(peak) then
              write(4,'(a9,2i4)') 'relocate ',ix,iy+1
              write(4,'("dot")')
            endif
            ipt1 = nx*(iy-1)+ix
            ipt2 = nx*iy+ix
            ipt3 = nx*iy+ix+1
            ipt4 = nx*(iy-1)+ix+1
            z1 = z(ipt1)
            z2 = z(ipt2)
            z3 = z(ipt3)
            z4 = z(ipt4)
            if(sqt) then
              z1 = sqrt(max(0.,z1))
              z2 = sqrt(max(0.,z2))
              z3 = sqrt(max(0.,z3))
              z4 = sqrt(max(0.,z4))
            endif
            if((z1*z2*z3*z4)==0.) cycle
            zhi = max(max(max(z1,z2),z3),z4)
            zlo = min(min(min(z1,z2),z3),z4)
            izlo = int((zlo-zmin)/dz+2.)-1
            izhi = int((zhi-zmin)/dz+1.)-1
            if(nz>0) then
              izlo = max(izlo,1)
              izhi = min(izhi,nz)
            else
              mzlo = 0
              do
                if(izlo>=nzabs) exit
                izlo = izlo/2
                mzlo = mzlo+1
              enddo
              izlo = izlo+mzlo*nzabs/2
              mzhi = 0
              do
                if(izhi>=nzabs) exit
                izhi = izhi/2
                mzhi = mzhi+1
              enddo
              izhi = izhi+mzhi*nzabs/2
              izlo = max(izlo,1)
            endif
            ncon = izhi-izlo+1
            if(ncon<1) cycle
            zc = (z1+z2+z3+z4)/4.
            if(transp) then
              xa = iy
              ya = ix
              dx = 0.
              dy = 1.
            else
              xa = ix
              ya = iy
              dx = 1.
              dy = 0.
            endif
            xc = xa+0.5
            yc = ya+0.5
            if(mirx) then
              xa = -xa
              xc = -xc
              dx = -dx
            endif
            if(miry) then
              ya = -ya
              yc = -yc
              dy = -dy
            endif
            ipta = ipt1
            za = z1
            jpt = 1
            do ia=1,2
              do ib=1,2
                iptb = ipta+jpt
                xb = xa+dx
                yb = ya+dy
                zb = z(iptb)
                if(sqt) zb = sqrt(max(0.,zb))
                zlo = min(min(za,zb),zc)
                zhi = max(max(za,zb),zc)
                izlo = int((zlo-zmin)/dz+2.)-1
                izhi = int((zhi-zmin)/dz+1.)-1
                if(nz>0) then
                  izlo = max(izlo,1)
                  izhi = min(izhi,nz)
                else
                  mzlo = 0
                  do
                    if(izlo>=nzabs) exit
                    izlo = izlo/2
                    mzlo = mzlo+1
                  enddo
                  izlo = izlo+mzlo*nzabs/2
                  mzhi = 0
                  do
                    if(izhi>=nzabs) exit
                    izhi = izhi/2
                    mzhi = mzhi+1
                  enddo
                  izhi = izhi+mzhi*nzabs/2
                  izlo = max(izlo,1)
                endif
                if(izhi>=izlo) then
                  do iz = izlo,izhi
                    start = .false.
                    if((nz>0).or.(iz<=nzabs)) then
                      zi = zmin+iz*dz
                    else
                      mz = (2*iz)/nzabs
                      zi = zmin+(iz-(mz-1)*nzabs/2)*(2**(mz-1))*dz
                    endif
                    if(za.ne.zb) then
                      dab = (zi-za)/(zb-za)
                      if((dab>=0.).and.(dab<1.)) then
                        start = .true.
                        x = xa+dab*(xb-xa)
                        y = ya+dab*(yb-ya)
                        write(4,'(a8,1p,2e12.3)') 'relocate',x,y
                      endif
                    endif
                    if(za.ne.zc) then
                      dac = (zi-za)/(zc-za)
                      if((dac>=0.).and.(dac<1.)) then
                        x = xa+dac*(xc-xa)
                        y = ya+dac*(yc-ya)
                        if(start) then
                          write(4,'(a4,1p,2e12.3)') 'draw',x,y
                          cycle
                        else
                          write(4,'(a8,1p,2e12.3)') 'relocate',x,y
                        endif
                      endif
                    endif
                    if(zb.ne.zc) then
                      dbc = (zi-zb)/(zc-zb)
                      if((dbc>=0.).and.(dbc<1.)) then
                        x = xb+dbc*(xc-xb)
                        y = yb+dbc*(yc-yb)
                        write(4,'(a4,1p,2e12.3)') 'draw',x,y
                      endif
                    endif
                  enddo
                endif
              enddo
              jpt = nx*jpt
              sdy = dy
              dy = xymir*dx
              dx = xymir*sdy
            enddo
            ipta = ipt3
            za = z(ipta)
            if(sqt) za = sqrt(max(0.,za))
            if(mirx) then
              xa = xa-1.
            else
              xa = xa+1.
            endif
            if(miry) then
              ya = ya-1.
            else
              ya = ya+1.
            endif
            dx = -dx
            dy = -dy
            jpt = -1
          enddo
        enddo
        write(4,'(a7)') 'angle 0'
      close (unit=4)

      end



      subroutine perspec(nx,ny,arr,scali,iwin)

        dimension arr(1)
        character*8 yn
        character*24 monfile
        logical ask,ok
        dimension yf(4096)

        if(iwin<1) return
        call pgslct(iwin)
        call pgsch(1.)
        call pgsubp(1,1)

        scale = scali
        ask = (scale==0.)
        if(scale==-1.) scale = 0.

        minx = 1
        maxx = nx
        miny = 1
        maxy = ny
        my = maxy-miny+1
        if(scale>0.) then
          zmn = 0.
          zmx = scale
        elseif(scale<0.) then
          zmn = scale
          zmx = -scale
        else
          zmx = 0.
          zmn = 0.
          do i = 1,nx*ny
            if(arr(i)<zmn) zmn = arr(i)
            if(arr(i)>zmx) zmx = arr(i)
          enddo
        endif
   10   if(ask) then
          if((zmx<(1.0e+06)).and.(zmx>1.0)) then
            print '(/" x,y,z range =",4i6,2f10.3,"  OK? (y)")', &
              minx,maxx,miny,maxy,zmn,zmx
          else
            print '(/" x,y,z range =",4i6,1p,2e10.2,"  OK? (y)")', &
              minx,maxx,miny,maxy,zmn,zmx
          endif
          read '(a8)', yn
          if(index(yn,'n')>0) then
            print '(" Enter x, y, z range")'
            read *, minx,maxx,miny,maxy,zmn,zmx
            my = maxy-miny+1
          endif
        endif

        xmn = minx
        xmx = maxx
        zrng = zmx-zmn
        scale = zrng/4.
        ymn = zmn/scale
        ymx = zmx/scale+my
        difmax = 0.01*my

        call pgenv(xmn,xmx,ymn,ymx,0,-2)

        do ix = 1,nx
          yf(ix) = ymn
        enddo
        do iy = 1,my
          ix = 1
          ipt = nx*(iy+miny-2)+ix
          x = ix
          ythis = iy
          y = ythis
          call pgmove(x,y)
          y = iy+arr(ipt)/scale
          call pgdraw(x,y)
          if(ythis>=yf(ix)) then
            ok = .true.
            yf(ix) = ythis
          else
            ok = .false.
          endif
          do ix = 2,nx
            ylast = ythis
            ipt = ipt+1
            if((ix==2).or.(ix==nx)) then
              amid = (arr(ipt-1)+arr(ipt))/2.
            else
!	cubic interpolation for mid point
              amid = (9.*(arr(ipt-1)+arr(ipt)) &
                -(arr(ipt-2)+arr(ipt+1)))/16.
            endif
            xmid = ix-0.5
            ymid = iy+amid/scale
            ythis = iy+arr(ipt)/scale
            if(ythis>=yf(ix)) then
!	this point is visible
              if(ok) then
!	and last point is visible
                if(abs((arr(ipt-1)+arr(ipt)) &
                  -(arr(ipt-2)+arr(ipt+1)))>difmax) &
                  call pgdraw(xmid,ymid)
                x = ix
                y = ythis
                call pgdraw(x,y)
              else
!	or last point is hidden
                ok = .true.
                x = ix-0.9
                y = 0.8*ylast+0.2*ymid
                call pgdraw(x,y)
                call pgmove(xmid,ymid)
                x = ix
                y = ythis
                call pgdraw(x,y)
              end if
              yf(ix) = ythis
            else
!	or this point is hidden
              if(ok) then
!	and last point is visible
                call pgdraw(xmid,ymid)
                ok = .false.
              else
!	or last point is hidden
                x = ix-0.9
                y = 0.8*ylast+0.2*ymid
                call pgdraw(x,y)
              endif
              x = ix-0.1
              y = 0.2*ymid+0.8*ythis
              call pgmove(x,y)
              x = ix
              y = ythis
              call pgdraw(x,y)
            end if
          enddo
        enddo

        if(.not.ask) return
        print '(" Replot screen plot? (n)")'
        read '(a8)', yn
        if(index(yn,'y')>0) go to 10

   39   print '(" Do you want a mongo plot file? (n)")'
        read '(a8)', yn
        if(index(yn,'y')==0) return
        print '(" Enter name of mongo file")'
        read '(a16)', monfile
        open(unit=4,file=monfile,status='new',access='sequential',err=39)
        write(4,*) 'reset'
        write(4,*) 'erase'
        write(4,*) 'expand 1.0'
        write(4,*) 'psfmode 2'
        write(4,'(a7,a16)') 'tlabel ',monfile
        write(4,'(a8,1p,2e12.3)') 'xlimits ',xmn,xmx
        write(4,'(a8,1p,2e12.3)') 'ylimits ',ymn,ymx
        write(4,*) 'window 1 1 1'
        write(4,*) 'expand 1.0'
        write(4,*) 'box 1 2 0 0'

        do ix = 1,nx
          yf(ix) = 0.
        enddo
        do iy = 1,my
          ix = 1
          ipt = nx*(iy+miny-2)+ix
          x = ix
          ythis = iy
          y = ythis
          write(4,'(a8,1p,2e12.3)') 'relocate',x,y
          y = iy+arr(ipt)/scale
          write(4,'(a4,1p,2e12.3)') 'draw',x,y
          if(ythis>=yf(ix)) then
            ok = .true.
            yf(ix) = ythis
          else
            ok = .false.
          endif
          do ix = 2,nx
            ylast = ythis
            ipt = ipt+1
            if((ix==2).or.(ix==nx)) then
              amid = (arr(ipt-1)+arr(ipt))/2.
            else
!	  cubic interpolation for mid point
              amid = (9.*(arr(ipt-1)+arr(ipt)) &
                     -(arr(ipt-2)+arr(ipt+1)))/16.
            endif
            xmid = ix-0.5
            ymid = iy+amid/scale
            ythis = iy+arr(ipt)/scale
            if(ythis>=yf(ix)) then
!	  this point is visible
              if(ok) then
!	  and last point is visible
                if(abs((arr(ipt-1)+arr(ipt)) &
                  -(arr(ipt-2)+arr(ipt+1)))>difmax) &
                  write(4,'(a4,1p,2e12.3)') 'draw',xmid,ymid
                  x = ix
                  y = ythis
                  write(4,'(a4,1p,2e12.3)') 'draw',x,y
                else
!	  or last point is hidden
                  ok = .true.
                  x = ix-0.9
                  y = 0.8*ylast+0.2*ymid
                  write(4,'(a4,1p,2e12.3)') 'draw',x,y
                  write(4,'(a8,1p,2e12.3)') 'relocate',xmid,ymid
                  x = ix
                  y = ythis
                  write(4,'(a4,1p,2e12.3)') 'draw',x,y
                end if
                yf(ix) = ythis
            else
!	  or this point is hidden
                if(ok) then
!	  and last point is visible
                  write(4,'(a4,1p,2e12.3)') 'draw',xmid,ymid
                  ok = .false.
                else
!	  or last point is hidden
                  x = ix-0.9
                  y = 0.8*ylast+0.2*ymid
                  write(4,'(a4,1p,2e12.3)') 'draw',x,y
                endif
                x = ix-0.1
                y = 0.2*ymid+0.8*ythis
                write(4,'(a8,1p,2e12.3)') 'relocate',x,y
                x = ix
                y = ythis
                write(4,'(a4,1p,2e12.3)') 'draw',x,y
            endif
          enddo
        enddo
        close (unit=4)

 90     return
        end


        subroutine histo(nx,ny,z,zpk,iwin)

          dimension z(1),h(1025),hsm(1025)

          if(iwin<1) return
          npt = nx*ny

          zmin = z(1)
          zmax = z(1)
          do ipt=1,npt
            zmin = amin1(zmin,z(ipt))
            zmax = amax1(zmax,z(ipt))
          enddo
          nz = min0(nx,1024)

          print '(" zmin, zmax, nz =",1p,2e10.2,i8)', &
            zmin,zmax,nz
!          read '(a8)', yn
!  24      if(index(yn,'n')>0) then
!            print '(" Enter zmin, dz, nz")'
!            read *, zmin,zmax,nz
!          endif
          dz = (zmax-zmin)/nz

          do iz = 1,nz
            h(iz) = 0.
          enddo
          i0 = 0

          do ipt = 1,npt
            if(z(ipt).ne.0.) then
                ibin = (z(ipt)-zmin)/dz+1.
                h(ibin) = h(ibin)+1.
            else
                i0 = i0+1
            endif
          enddo
          ibin0 = -zmin/dz+1.
          print '(" histogram(",i3," = 0) =",i5)', ibin0,i0

          ng = 1
          s = 0.
          call graph(nz,ng,h,s,iwin)

!          print '(" Replot histogram? (n)")'
!          read '(a8)', yn
!          if(index(yn,'y')>0) go to 20

          nsm = 64/(npt/nz)
          nsm = max(1,(min(8,nsm)))
          sm = 0.
          do ibin = 1,nz+2*nsm
            if(ibin<=nz) sm = sm+h(ibin)
            jbin = ibin-nsm
            kbin = jbin-nsm
            if((jbin>=1).and.(jbin<=nz)) hsm(jbin) = sm
            if(kbin>=1) sm = sm-h(kbin)
          enddo
          sm = 0.
          smx = 0.
          do ibin = 1,nz+2*nsm
            if(ibin<=nz) sm = sm+hsm(ibin)
            jbin = ibin-nsm
            kbin = jbin-nsm
            if(sm>smx) then
              smx = sm
              jmx = jbin
            endif
            if(kbin>=1) sm = sm-hsm(kbin)
          enddo
          zpk = zmin+dz*(jmx-0.5)
          print '(" Peak of smoothed histogram at ",1p,d10.2)', zpk

          return
        end




        subroutine graph(nx,ny,array,scale,iwin)

        character*8 yn
        character*24 monfile
        dimension array(1)

        if(iwin<1) return
        call pgslct(iwin)
        call pgsch(1.)
        call pgsubp(1,1)

 10     xmn = 1.
        xmx = nx
        if(scale>0.) then
          ymn = 0.
          ymx = scale
        elseif(scale<0.) then
          ymn = scale
          ymx = -scale
        else
          ymx = 0.
          ymn = 0.
          do  i=1,nx*ny
            if(array(i)<ymn) ymn = array(i)
            if(array(i)>ymx) ymx = array(i)
          enddo
        endif
!        if((xmx<(1.0e+06)).and.(ymx<(1.0e+06))) then
!          print'(/" xmn,xmx,ymn,ymx=",4f10.3,"  OK? (y)")', &
!            xmn,xmx,ymn,ymx
!        else
!          print'(/" xmn,xmx,ymn,ymx=",1p,4e10.2,"  OK? (y)")', &
!            xmn,xmx,ymn,ymx
!        endif
!        read '(a8)', yn
!        if(index(yn,'n')>0) then
!          print '(" Enter xmn, xmx, ymn, ymx")'
!          read *, xmn,xmx,ymn,ymx
!        endif
        xrng = xmx-xmn
        yrng = ymx-ymn

        call pgenv(xmn,xmx,ymn,ymx,0,0)

        ipt = 0
        do iy = 1,ny
          ltype = mod(iy-1,4)+1
          call pgsls(ltype)
          x = 1
          y = array(ipt+1)
          y0 = y
          call pgmove(x,y)
          do ix = 1,nx
            ipt = ipt+1
            x = ix
            if((x<xmn).or.(x>xmx)) cycle
            y = array(ipt)
            call pgdraw(x,y)
            y0 = y
          enddo
        enddo

!        print '(" Replot plot? (n)")'
!        read '(a8)', yn
!        if(index(yn,'y')>0) go to 10

   39   continue
!        print '(" Do you want a mongo plot input file? (n)")'
!        read '(a8)', yn
        yn = 'n'
        if(index(yn,'y')==0) return
        print '(" Enter name of mongo file")'
        read '(a16)', monfile
        open(unit=4,file=monfile,status='new',access='sequential',err=39)
        write(4,'(a5)') 'reset'
        write(4,'(a5)') 'erase'
        write(4,'(a10)') 'expand 1.0'
        write(4,'(a9)') 'psfmode 2'
        write(4,'(a8,1p,2e12.3)') 'xlimits ',xmn,xmx
        write(4,'(a8,1p,2e12.3)') 'ylimits ',ymn,ymx
        write(4,'(a12)') 'window 1 1 1'
        write(4,'(a10)') 'expand 1.0'
        write(4,'(a11)') 'box 1 2 0 0'
        xrng = xmx-xmn
        yrng = ymx-ymn
        write(4,'(a7,1p,2e12.3)') 'xlimits',xmn/xrng,xmx/xrng
        write(4,'(a7,1p,2e12.3)') 'ylimits',ymn/yrng,ymx/yrng

        ipt = 0
        do iy = 1,ny
          ltype = mod(iy-1,6)
          write(4,'(a5,i2)') 'ltype',ltype
          x = 1
          y = array(ipt+1)
          y0 = y
          write(4,'(a8,1p,2e12.3)') 'relocate',x/xrng,y/yrng
          do ix = 1,nx
            ipt = ipt+1
            x = ix
            if((x<xmn).or.(x>xmx)) cycle
            y = array(ipt)
            write(4,'(a4,1p,2e12.3)') 'draw',x/xrng,y/yrng
            y0 = y
          enddo
        enddo

        close(unit=4)

        return
        end




        subroutine scatter(npt,xarr,yarr,scale,datafile,iwin)

        character*8 yn
        dimension xarr(1),yarr(1)
        character*24 monfile,datafile
        logical mongo

        length = index(datafile,' ')-1
        mongo = (length>0)

        if(iwin<1) return
        call pgslct(iwin)
        call pgsch(1.)
        call pgsubp(1,1)

        if(scale>0.) then
          xmn = 0.
          xmx = scale
          ymn = 0.
          ymx = scale
        elseif(scale<0.) then
          xmn = scale
          xmx = -scale
          ymn = scale
          ymx = -scale
        else
          xmn = 0.
          xmx = 0.
          ymn = 0.
          ymx = 0.
          do i=1,npt
            if(xarr(i)<xmn) xmn = xarr(i)
            if(xarr(i)>xmx) xmx = xarr(i)
            if(yarr(i)<ymn) ymn = yarr(i)
            if(yarr(i)>ymx) ymx = yarr(i)
          enddo
        endif
   12   if((xmx<(1.0e+06)).and.(ymx<(1.0e+06))) then
          print'(/" xmn,xmx,ymn,ymx=",4f10.3,"  OK? (y)")', &
            xmn,xmx,ymn,ymx
        else
          print'(/" xmn,xmx,ymn,ymx=",1p,4e10.2,"  OK? (y)")', &
            xmn,xmx,ymn,ymx
        endif
        read '(a8)', yn
        if(index(yn,'n')>0) then
          print '(" Enter xmn, xmx, ymn, ymx")'
          read *, xmn,xmx,ymn,ymx
        endif
        dxpt = (xmx-xmn)/1000.

        call pgenv(xmn,xmx,ymn,ymx,0,0)

        do ipt = 1,npt
          x = xarr(ipt)
          y = yarr(ipt)
          if((x>=xmn).and.(x<=xmx) &
            .and.(y>=ymn).and.(y<=ymx)) then
            call pgmove(x-dxpt,y)
            call pgdraw(x+dxpt,y)
            call pgmove(x,y-dxpt)
            call pgdraw(x,y+dxpt)
          endif
        enddo

        print '(" Replot plot? (n)")'
        read '(a8)', yn
        if(index(yn,'y')>0) go to 12

        if(mongo) then

          monfile = datafile(:length)//'.inp'
          open(unit=4,file=monfile,access='sequential',err=39)
          monfile = datafile(1:length)//'.mon'
          write(4,'(a5,a24)') 'data ',monfile
          write(4,'(a5)') 'reset'
          write(4,'(a5)') 'erase'
          write(4,'(a9)') 'psfmode 2'
          write(4,'(a7,a)') 'tlabel ',datafile(1:length)

          ni = maxx-minx+1
          write(4,'(a8,1p,2e12.3)') 'ylimits ',ymn,ymx
          write(4,'("window 1 1 1")')
          write(4,'(a8,2es12.3)') 'xlimits ',xmn,xmx
          write(4,'(a11)') 'box 1 2 0 0'
          write(4,'("lines ",2i8)') 1,npt
          write(4,'("ycolumn 2")')
          write(4,'("xcolumn 1")')
          write(4,'("expand 0.5")')
          write(4,'("ptype 4 1")')
          write(4,'("points")')
          close(unit=4)

          monfile = datafile(:length)//'.mon'
          open(unit=4,file=monfile,access='sequential',err=39)
          do ipt = 1,npt
            write(4,'(2es12.3)') xarr(ipt),yarr(ipt)
          enddo
          close(unit=4)

        endif

   39   return
        end



        subroutine plots(mx,nx,my,arr,scale,datafile,iwin)
!	mx: x dim of arr(mx,5), which is interpreted as arr(nx,ny,5)
!	ny = |my|, with my < 0 to suppress xlabels
!	iz = 1: x axis, iz = 2,3: y axis, iz = 4: atmo plot

        character*8 yn
        character*24 monfile,datafile
        logical ask,save,nolabel
        dimension arr(mx,5)
        real, parameter :: big = 1.0e+31

        if(iwin<1) return
        call pgslct(iwin)
        ask = (scale==0.)
        if(scale==-1.) scale = 0.
        istart = index(datafile,'/')+1
        length = index(datafile,' ')-1
        save = (length>0)
        nolabel = (my<0)
        ny = iabs(my)
        nv = int(sqrt(float(ny)))
        nu = nv*nx
        nw = (ny-1)/nv+1

 10     minx = 1
        maxx = nx
        if(scale>0.) then
          ymn = 0.
          ymx = scale
        elseif(scale<0.) then
          ymn = scale
          ymx = -scale
        else
          ymx = 0.
          ymn = 0.
          do iy=1,ny
            do ix=1,nx-2
              i = nx*(iy-1)+ix
              if((arr(i,2)<ymn).and.(arr(i+1,2)<ymn) &
                 .and.(arr(i+2,2)<ymn)) &
                ymn = max(arr(i,2),arr(i+1,2),arr(i+2,2))
              if((arr(i,2)>ymx).and.(arr(i+1,2)>ymx) &
                 .and.(arr(i+2,2)>ymx)) &
                ymx = 1.02*min(arr(i,2),arr(i+1,2),arr(i+2,2))
            enddo
          enddo
          scale = ymx
          if(ymn<(-ymx/10.)) scale = -scale
        endif
        if(ask.or.(ymn>=ymx).or.(ymn.ne.ymn).or.(ymx.ne.ymx) &
           .or.(ymn<(-big)).or.(ymx>big)) then
          if(ymx<(1.0e+06)) then
            print '(/" xmn,xmx,ymn,ymx=",2i6,2f10.3,"  OK? (y)")', &
              minx,maxx,ymn,ymx
          else
            print '(/" xmn,xmx,ymn,ymx=",2i6,1p,2e10.2,"  OK? (y)")', &
              minx,maxx,ymn,ymx
          endif
          read '(a8)', yn
          if(index(yn,'n')>0) then
            print '(" Enter xmn, xmx, ymn, ymx")'
            read *, minx,maxx,ymn,ymx
          endif
        endif
        midx = (minx+maxx)/2

        call pgsubp(nv,nw)
        if(nolabel) then
          call pgsch(0.)
        else
          call pgsch(1.)
        endif

        do iy = 1,ny
          ix0 = nx*(iy-1)
          xmn = arr(ix0+minx,1)
          xmx = arr(ix0+maxx,1)
          if(xmn==xmx) then
            do ix = minx,maxx
              arr(ix0+ix,1) = ix
            enddo
            xmn = arr(ix0+minx,1)
            xmx = arr(ix0+maxx,1)
          endif
          if(xmn>xmx) then
            xxx = xmn
            xmn = xmx
            xmx = xxx
          endif
          dx2 = (arr(ix0+midx+1,1)-arr(ix0+midx,1))/2.
          if(nolabel) then
            call pgenv(xmn,xmx,ymn,ymx,0,-1)
          else
            call pgenv(xmn,xmx,ymn,ymx,0,0)
          endif
          x = arr(ix0+minx,1)
          y = arr(ix0+minx,2)
          if(.not.(0.*y==0.)) y = 0.
          call pgmove(x,y)
          do ix = minx+1,maxx
            x = arr(ix0+ix,1)-dx2
            if((x<xmn).or.(x>xmx)) cycle
            call pgdraw(x,y)
            y = arr(ix0+ix,2)
            if(.not.(0.*y==0.)) y = 0.
            call pgdraw(x,y)
          enddo
          x = arr(ix0+maxx,1)
          if((x>=xmn).and.(x<=xmx)) call pgdraw(x,y)
          x = arr(ix0+minx,1)
          y = arr(ix0+minx,3)
          if(.not.(0.*y==0.)) y = 0.
          call pgmove(x,y)
          do ix = minx+1,maxx
            x = arr(ix0+ix,1)-dx2
            if((x<xmn).or.(x>xmx)) cycle
            call pgdraw(x,y)
            y = arr(ix0+ix,3)
            if(.not.(0.*y==0.)) y = 0.
            call pgdraw(x,y)
          enddo
          x = arr(ix0+maxx,1)
          if((x>=xmn).and.(x<=xmx)) call pgdraw(x,y)
          x = arr(ix0+maxx,1)
          if((x>=xmn).and.(x<=xmx)) call pgdraw(x,y)
          x = arr(ix0+minx,1)
          y = arr(ix0+minx,4)*(ymx/2.)
          if(.not.(0.*y==0.)) y = 0.
          call pgmove(x,y)
          do ix = minx+1,maxx
            x = arr(ix0+ix,1)-dx2
            if((x<xmn).or.(x>xmx)) cycle
            call pgdraw(x,y)
            y = arr(ix0+ix,4)*(ymx/2.)
            if(.not.(0.*y==0.)) y = 0.
            call pgdraw(x,y)
          enddo
          x = arr(ix0+maxx,1)
          if((x>=xmn).and.(x<=xmx)) call pgdraw(x,y)
        enddo

        if(ask) then
          print '(" Replot plot? (n)")'
          read '(a8)', yn
          if(index(yn,'y')>0) go to 10
        endif

        if(save) then
          monfile = datafile(:length)//'.inp'
          open(unit=4,file=monfile,access='sequential',err=39)
          monfile = datafile(istart:length)//'.mon'
          write(4,'(a5,a24)') 'data ',monfile
          write(4,'(a5)') 'reset'
          write(4,'(a5)') 'erase'
          write(4,'(a10)') 'expand 0.8'
          write(4,'(a)') 'psfmode 2'
          write(4,'(a7,a)') 'tlabel ',datafile(istart:length)

          ni = maxx-minx+1
          write(4,'(a18)') 'submargins 0.0 0.5'
          write(4,'(a8,1p,2e12.3)') 'ylimits ',ymn,ymx
          do iw = 1,nw
            xmn = arr((iw-1)*nu+1,1)
            if((iw*nv)<=ny) then
              xmx = arr(iw*nu,1)
            else
              xmx = xmx+xmn-arr((iw-2)*nu+1,1)
            endif
            write(4,'(a6,3i3)') 'window', 1,nw,(nw-iw+1)
            write(4,'(a8,2f12.6)') 'xlimits ',xmn,xmx
            write(4,'(a10)') 'expand 0.5'
            write(4,'(a11)') 'box 1 2 0 0'
            do iv = 1,nv
              iy = (iw-1)*nv+iv
              if(iy>ny) cycle
              ix2 = iy*nx
              ix1 = ix2-nx+1
              write(4,'("lines ",2i8)') ix1,ix2
              write(4,'("ycolumn 2")')
              write(4,'("xcolumn 1")')
              write(4,'("histogram")')
            enddo
          enddo
          close(unit=4)

          monfile = datafile(:length)//'.mon'
          open(unit=4,file=monfile,access='sequential',err=39)
          do iy = 1,ny
            ix0 = nx*(iy-1)
            do ix = 1,nx
              x = arr(ix0+ix,1)
              y = arr(ix0+ix,2)
              z = arr(ix0+ix,3)
              a = arr(ix0+ix,4)
              if(.false.) then
                w = 1.0e+04/x
                write(4,'(f12.4,3es12.3,f12.6)') x,y,z,a,w
              else
                write(4,'(f12.4,4es12.3)') x,y,z,a,arr(ix0+ix,5)
              endif
            enddo
          enddo
          close(unit=4)

        endif

   39   return
        end




        subroutine plot4(mx,nx,ny,arr,scale,datafile,iwin)
!	mx: x dim of arr(mx,5), which is interpreted as arr(nx,ny,5)
!	nx = points in each plot window, ny = number of plot windows
!	iz = 1: x axis, iz = 2-5: y axis for 4 plots

        character*8 yn
        character*24 monfile,datafile
        logical ask,save,nolabel
        dimension arr(mx,5)
        real, parameter :: big = 1.0e+31

        if(iwin<1) return
        call pgslct(iwin)
        ask = (scale==0.)
        if(scale==-1.) scale = 0.
        istart = index(datafile,'/')+1
        length = index(datafile,' ')-1
        save = (length>0)
        nv = int(sqrt(float(ny)))
        nu = nv*nx
        nw = (ny-1)/nv+1

 10     minx = 1
        maxx = nx
        if(scale>0.) then
          ymn = 0.
          ymx = scale
        elseif(scale<0.) then
          ymn = scale
          ymx = -scale
        else
          ymx = 0.
          ymn = 0.
          do iy=1,ny
            do ix=1,nx-2
              i = nx*(iy-1)+ix
              if((arr(i,2)<ymn).and.(arr(i+1,2)<ymn) &
                 .and.(arr(i+2,2)<ymn)) &
                ymn = max(arr(i,2),arr(i+1,2),arr(i+2,2))
              if((arr(i,2)>ymx).and.(arr(i+1,2)>ymx) &
                 .and.(arr(i+2,2)>ymx)) &
                ymx = 1.02*min(arr(i,2),arr(i+1,2),arr(i+2,2))
            enddo
          enddo
          scale = ymx
          if(ymn<(-ymx/10.)) scale = -scale
        endif
        if(ask.or.(ymn>=ymx).or.(ymn.ne.ymn).or.(ymx.ne.ymx) &
           .or.(ymn<(-big)).or.(ymx>big)) then
          if(ymx<(1.0e+06)) then
            print '(/" xmn,xmx,ymn,ymx=",2i6,2f10.3,"  OK? (y)")', &
              minx,maxx,ymn,ymx
          else
            print '(/" xmn,xmx,ymn,ymx=",2i6,1p,2e10.2,"  OK? (y)")', &
              minx,maxx,ymn,ymx
          endif
          read '(a8)', yn
          if(index(yn,'n')>0) then
            print '(" Enter xmn, xmx, ymn, ymx")'
            read *, minx,maxx,ymn,ymx
          endif
        endif
        midx = (minx+maxx)/2

        call pgsubp(nv,nw)
        if(nolabel) then
          call pgsch(0.)
        else
          call pgsch(1.)
        endif

        do iy = 1,ny
          ix0 = nx*(iy-1)
          xmn = arr(ix0+minx,1)
          xmx = arr(ix0+maxx,1)
          if(xmn==xmx) then
            do ix = minx,maxx
              arr(ix0+ix,1) = ix
            enddo
            xmn = arr(ix0+minx,1)
            xmx = arr(ix0+maxx,1)
          endif
          if(xmn>xmx) then
            xxx = xmn
            xmn = xmx
            xmx = xxx
          endif
          dx2 = (arr(ix0+midx+1,1)-arr(ix0+midx,1))/2.
          call pgenv(xmn,xmx,ymn,ymx,0,0)
          x = arr(ix0+minx,1)
          y = arr(ix0+minx,2)
          if(.not.(0.*y==0.)) y = 0.
          call pgsci(1)
          call pgmove(x,y)
          do ix = minx+1,maxx
            x = arr(ix0+ix,1)-dx2
            if((x<xmn).or.(x>xmx)) cycle
            call pgdraw(x,y)
            y = arr(ix0+ix,2)
            if(.not.(0.*y==0.)) y = 0.
            call pgdraw(x,y)
          enddo
          call pgsci(3)
          x = arr(ix0+maxx,1)
          if((x>=xmn).and.(x<=xmx)) call pgdraw(x,y)
          x = arr(ix0+minx,1)
          y = arr(ix0+minx,3)
          if(.not.(0.*y==0.)) y = 0.
          call pgmove(x,y)
          do ix = minx+1,maxx
            x = arr(ix0+ix,1)-dx2
            if((x<xmn).or.(x>xmx)) cycle
            call pgdraw(x,y)
            y = arr(ix0+ix,3)
            if(.not.(0.*y==0.)) y = 0.
            call pgdraw(x,y)
          enddo
          x = arr(ix0+maxx,1)
          if((x>=xmn).and.(x<=xmx)) call pgdraw(x,y)
          x = arr(ix0+maxx,1)
          if((x>=xmn).and.(x<=xmx)) call pgdraw(x,y)
          x = arr(ix0+minx,1)
          y = arr(ix0+minx,4)
          if(.not.(0.*y==0.)) y = 0.
          call pgmove(x,y)
          do ix = minx+1,maxx
            x = arr(ix0+ix,1)-dx2
            if((x<xmn).or.(x>xmx)) cycle
            call pgdraw(x,y)
            y = arr(ix0+ix,4)
            if(.not.(0.*y==0.)) y = 0.
            call pgdraw(x,y)
          enddo
          x = arr(ix0+maxx,1)
          if((x>=xmn).and.(x<=xmx)) call pgdraw(x,y)
          x = arr(ix0+maxx,1)
          if((x>=xmn).and.(x<=xmx)) call pgdraw(x,y)
          x = arr(ix0+minx,1)
          y = arr(ix0+minx,5)
          if(.not.(0.*y==0.)) y = 0.
          call pgmove(x,y)
          do ix = minx+1,maxx
            x = arr(ix0+ix,1)-dx2
            if((x<xmn).or.(x>xmx)) cycle
            call pgdraw(x,y)
            y = arr(ix0+ix,5)
            if(.not.(0.*y==0.)) y = 0.
            call pgdraw(x,y)
          enddo
          x = arr(ix0+maxx,1)
          if((x>=xmn).and.(x<=xmx)) call pgdraw(x,y)
          call pgsci(1)
        enddo

        if(ask) then
          print '(" Replot plot? (n)")'
          read '(a8)', yn
          if(index(yn,'y')>0) go to 10
        endif

        if(save) then
          monfile = datafile(:length)//'.inp'
          open(unit=4,file=monfile,access='sequential',err=39)
          monfile = datafile(istart:length)//'.mon'
          write(4,'(a5,a24)') 'data ',monfile
          write(4,'(a5)') 'reset'
          write(4,'(a5)') 'erase'
          write(4,'(a10)') 'expand 0.8'
          write(4,'(a)') 'psfmode 2'
          write(4,'(a7,a)') 'tlabel ',datafile(istart:length)

          ni = maxx-minx+1
          write(4,'(a18)') 'submargins 0.0 0.5'
          write(4,'(a8,1p,2e12.3)') 'ylimits ',ymn,ymx
          do iw = 1,nw
            xmn = arr((iw-1)*nu+1,1)
            if((iw*nv)<=ny) then
              xmx = arr(iw*nu,1)
            else
              xmx = xmx+xmn-arr((iw-2)*nu+1,1)
            endif
            write(4,'(a6,3i3)') 'window', 1,nw,(nw-iw+1)
            write(4,'(a8,2f12.6)') 'xlimits ',xmn,xmx
            write(4,'(a10)') 'expand 0.5'
            write(4,'(a11)') 'box 1 2 0 0'
            do iv = 1,nv
              iy = (iw-1)*nv+iv
              if(iy>ny) cycle
              ix2 = iy*nx
              ix1 = ix2-nx+1
              write(4,'("lines ",2i8)') ix1,ix2
              write(4,'("ycolumn 2")')
              write(4,'("xcolumn 1")')
              write(4,'("histogram")')
            enddo
          enddo
          close(unit=4)

          monfile = datafile(:length)//'.mon'
          open(unit=4,file=monfile,access='sequential',err=39)
          do iy = 1,ny
            ix0 = nx*(iy-1)
            do ix = 1,nx
              x = arr(ix0+ix,1)
              y = arr(ix0+ix,2)
              z = arr(ix0+ix,3)
              a = arr(ix0+ix,4)
              b = arr(ix0+ix,5)
              if(.false.) then
                w = 1.0e+04/x
                write(4,'(f12.4,3es12.3,f12.6)') x,y,z,a,b,w
              else
                write(4,'(f12.4,4es12.3)') x,y,z,a,b
              endif
            enddo
          enddo
          close(unit=4)

        endif

   39   return
        end




        subroutine plotg(mx,nx,ny,nz,arr,scale,iwin)
!	2-dimensional grid of plots
!	mx: x dim of arr(mx,4), which is interpreted as arr(nx,ny,nz,4)
!	i4 = 1: x axis, i4 = 2: y axis, i4 = 3,4: store, but don''t plot

        character*8 yn
        character*24 monfile,datafile
        logical ask,axes
        dimension arr(mx,4)
        real, parameter :: big = 1.0e+31

        if(iwin<1) return
        call pgslct(iwin)
        call pgsch(0.)
        ask = (scale==0.)
        if(scale==-1.) scale = 0.
        axes = (nx>0)
        nx = iabs(nx)

 10     minx = 1
        maxx = nx
        if(scale>0.) then
          ymn = 0.
          ymx = scale
        elseif(scale<0.) then
          ymn = scale
          ymx = -scale
        else
          ymx = 0.
          ymn = 0.
          do iz=1,nz
            do iy=1,ny
              i = nx*ny*(iz-1)+nx*(iy-1)
              do ix=1,nx-2
                i = i+1
                if((arr(i,2)<ymn).and.(arr(i+1,2)<ymn) &
                   .and.(arr(i+2,2)<ymn)) &
                  ymn = max(arr(i,2),arr(i+1,2),arr(i+2,2))
                if((arr(i,2)>ymx).and.(arr(i+1,2)>ymx) &
                    .and.(arr(i+2,2)>ymx)) &
                  ymx = 1.02*min(arr(i,2),arr(i+1,2),arr(i+2,2))
              enddo
            enddo
          enddo
          scale = ymx
          if(ymn<(-ymx/10.)) scale = -scale
        endif
        if(ask.or.(ymn.ne.ymn).or.(ymx.ne.ymx) &
           .or.(ymn<(-big)).or.(ymx>big)) then
          if(ymx<(1.0e+06)) then
            print '(/" xmn,xmx,ymn,ymx=",2i6,2f10.3,"  OK? (y)")', &
              minx,maxx,ymn,ymx
          else
            print '(/" xmn,xmx,ymn,ymx=",2i6,1p,2e10.2,"  OK? (y)")', &
              minx,maxx,ymn,ymx
          endif
          read '(a8)', yn
          if(index(yn,'n')>0) then
            print '(" Enter xmn, xmx, ymn, ymx")'
            read *, minx,maxx,ymn,ymx
          endif
        endif

        call pgsubp(ny,nz)

        do iz = nz,1,-1
          do iy = 1,ny
            ix0 = nx*ny*(iz-1)+nx*(iy-1)
            xmn = arr(ix0+minx,1)
            xmx = arr(ix0+maxx,1)
            if(xmx==xmn) print '(3i8,1p,2e10.3)', ix0,minx,maxx,xmn,xmx
            if(xmn>xmx) then
              xxx = xmn
              xmn = xmx
              xmx = xxx
            endif
            if(axes) then
              call pgenv(xmn,xmx,ymn,ymx,0,-1)
            else
              call pgenv(xmn,xmx,ymn,ymx,0,-2)
            endif
            x = arr(ix0+minx,1)
            y = arr(ix0+minx,2)
            call pgmove(x,y)
            do ix = minx,maxx
              x = arr(ix0+ix,1)
              if((x<xmn).or.(x>xmx)) cycle
              call pgdraw(x,y)
              y = arr(ix0+ix,2)
              call pgdraw(x,y)
            enddo
          enddo
        enddo

        if(ask) then
          print '(" Replot plot? (n)")'
          read '(a8)', yn
          if(index(yn,'y')>0) go to 10

          print '(" Do you want a mongo plot file? (n)")'
          read '(a8)', yn
          if(index(yn,'y')==0) return
          print '(" Enter name of mongo file")'
          read '(a16)', datafile
          istart = index(datafile,'/')+1
          length = index(datafile,' ')-1

          monfile = datafile(:length)//'.inp'
          open(unit=4,file=monfile,access='sequential',err=39)
          monfile = datafile(istart:length)//'.mon'
          write(4,'(a5,a24)') 'data ',monfile
          write(4,'(a5)') 'reset'
          write(4,'(a5)') 'erase'
          write(4,'(a10)') 'expand 0.8'
          write(4,'(a)') 'psfmode 2'
          write(4,'(a7,a)') 'tlabel ',datafile(istart:length)

          ni = maxx-minx+1
          write(4,'(a18)') 'submargins 0.0 0.0'
          write(4,'(a8,1p,2e12.3)') 'ylimits ',ymn,ymx
          do iz = 1,nz
            do iy = 1,ny
              ipl = ny*(iz-1)+iy
              ix0 = nx*ny*(iz-1)+nx*(iy-1)
              ix1 = ix0+1
              ix2 = ix0+nx
              xmn = arr(ix1,1)
              xmx = arr(ix2,1)
              write(4,'(a6,3i4)') 'window', ny,nz,iy+ny*(iz-1)
              write(4,'(a8,2f12.6)') 'xlimits ',xmn,xmx
              write(4,'(a10)') 'expand 0.5'
              if(axes) then
                if(iy==1) then
                  if(iz==1) then
                    write(4,'(a11)') 'box 1 2 5 5'
                  else
                    write(4,'(a11)') 'box 0 2 5 5'
                  endif
                elseif(iz==1) then
                  write(4,'(a11)') 'box 1 0 5 5'
                else
                  write(4,'(a11)') 'box 0 0 5 5'
                endif
              endif
              write(4,'("lines ",2i8)') ix1,ix2
              write(4,'("ycolumn 2")')
              write(4,'("xcolumn 1")')
              write(4,'("histogram")')
            enddo
          enddo
          close(unit=4)

          monfile = datafile(:length)//'.mon'
          open(unit=4,file=monfile,access='sequential',err=39)
          do iz = 1,nz
            do iy = 1,ny
              ix0 = nx*ny*(iz-1)+nx*(iy-1)
              do ix = 1,nx
                ixi = ix0+ix
                x = arr(ixi,1)
                y = arr(ixi,2)
                z = arr(ixi,3)
                a = arr(ixi,4)
                if(y==0.) z = 0.
                write(4,'(f12.3,1p,3e12.3)') x,y,z,a
              enddo
            enddo
          enddo
          close(unit=4)
        endif

        return

   39   print '(" Error opening mongo files")'
        return

      end




      subroutine getcurs(x,y,ikey,ierr)

        integer pgcurs,pgband
        character*1 c

        x = 0.
        y = 0.
!        ierr = pgcurs(x,y,c)
        ierr = pgband(7,1,x,y,x,y,c)
        ierr = ierr-1
        if(ierr==0) then
          read(unit=c,fmt='(i1)',err=90) ikey
        else
          print '(" Error returned by pgcurs")'
          ikey = -1
        endif
        return

 90     ikey = -1
        return

      end



      function plmin(pl,nx,x,nsmoo)

        real pl(nx,2)
        real plsm(16384)
        character*8 yn

   10   nsmoo = min(nsmoo,1)
        xnsm = nsmoo
        nsm1 = (nsmoo+1)/2
        nsm2 = nsmoo/2
        sum = nsm1*pl(1,2)
        do ix = 1,nsm2
          sum = sum+pl(ix,2)
        enddo
        do ix = 1,nx
          if(nsmoo<=1) then
            sum = xnsm*pl(ix,2)
          elseif(ix+nsm2>nx) then
            sum = sum+pl(nx,2)-pl(ix-nsm1,2)
          elseif(ix-nsm1<1) then
            sum = sum+pl(ix+nsm2,2)-pl(1,2)
          else
            sum = sum+pl(ix+nsm2,2)-pl(ix-nsm1,2)
          endif
          plsm(ix) = sum/xnsm
        enddo

        do ix = 1,nx
          if(pl(ix,1)>x) then
            i0 = ix
            exit
          endif
        enddo

        if((i0<=1).or.(i0>=nx)) then
          print '(" x out of range in plmin",2i6,1p,2e10.2)', &
            i0,nx,x,pl(nx,1)
          plmin = x
          return
        endif

        if(plsm(i0-1)<plsm(i0)) i0 = i0-1
        if(plsm(i0-1)<plsm(i0)) i0 = i0-1
        if(plsm(i0+1)<plsm(i0)) i0 = i0+1
        i1 = i0-1
        i2 = i0
        i3 = i0+1
        pl1 = plsm(i1)
        pl2 = plsm(i2)
        pl3 = plsm(i3)
        if((pl1<pl2).or.(pl3<pl2)) then
          print '(" Point chosen not at a minimum.", &
            "  Change smoothing? ",$)'
          read '(a8)', yn
          if(index(yn,'y')>0) then
            print '(" Enter nsmooth ")'
            read *, nsmoo
            goto 10
          else
            print '(" Find maximum? ",$)'
            read '(a8)', yn
            if(index(yn,'y')>0) then
              plmin = plmax(pl,nx,x,nsmoo)
              return
            endif
            print '(" Taking cursor position")'
            plmin = x
            return
          endif
        endif

        di = 0.5*(pl1-pl3)/(pl1-2.*pl2+pl3)
        pl1 = pl(i1,1)
        pl2 = pl(i2,1)
        pl3 = pl(i3,1)
        if(mod(nsmoo,2)==0) di = di+0.5
        plmin = pl2+di*(pl3-pl1)/2.
        print '(" Minimum at pixel ",f6.2," = ",f8.3)', i0+di,plmin

        return
      end


      function plmax(pl,nx,x,nsmoo)

        real pl(nx,2)
        real plsm(16384)
        character*8 yn

   10   nsm1 = (nsmoo+1)/2
        nsm2 = nsmoo/2
        sum = nsm1*pl(1,2)
        do ix = 1,nsm2
          sum = sum+pl(ix,2)
        enddo
        do ix = 1,nx
          if(nsmoo<=1) then
            sum = nsmoo*pl(ix,2)
          elseif(ix+nsm2>nx) then
            sum = sum+pl(nx,2)-pl(ix-nsm1,2)
          elseif(ix-nsm1<1) then
            sum = sum+pl(ix+nsm2,2)-pl(1,2)
          else
            sum = sum+pl(ix+nsm2,2)-pl(ix-nsm1,2)
          endif
          plsm(ix) = sum/nsmoo
        enddo

        do ix = 1,nx
          if(pl(ix,1)>x) then
            i0 = ix
            exit
          endif
        enddo

        if((i0<=1).or.(i0>=nx)) then
          print '(" x out of range in plmax",2i6,1p,2e10.2)', &
            i0,nx,x,pl(nx,1)
          plmax = x
          return
        endif

        if(plsm(i0-1)>plsm(i0)) i0 = i0-1
        if(plsm(i0-1)>plsm(i0)) i0 = i0-1
        if(plsm(i0+1)>plsm(i0)) i0 = i0+1
        i1 = i0-1
        i2 = i0
        i3 = i0+1
        pl1 = plsm(i1)
        pl2 = plsm(i2)
        pl3 = plsm(i3)
        if((pl1>pl2).or.(pl3>pl2)) then
          print '(" Point chosen not at a maximum.", &
            "  Change smoothing? ",$)'
          read '(a8)', yn
          if(index(yn,'y')>0) then
            print '(" Enter nsmooth ")'
            read *, nsmoo
            goto 10
          else
            print '(" Taking cursor position")'
            plmax = x
            return
          endif
        endif

        di = 0.5*(pl1-pl3)/(pl1-2.*pl2+pl3)
        pl1 = pl(i1,1)
        pl2 = pl(i2,1)
        pl3 = pl(i3,1)
        if(mod(nsmoo,2)==0) di = di+0.5
        plmax = pl2+di*(pl3-pl1)/2.
        print '(" Maximum at pixel ",f6.2," = ",f8.3)', i0+di,plmax

        return
      end

