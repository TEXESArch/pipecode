      module ius
        integer, parameter :: iuwno = 10
        integer, parameter :: iups = 11
        integer, parameter :: iupl = 12
        integer, parameter :: iurawh = 13
        integer, parameter :: iuredh = 14
        integer, parameter :: iusumh = 15
        integer, parameter :: iurawd = 17
        integer, parameter :: iuredd = 18
        integer, parameter :: iuatmo = 20
        integer, parameter :: iufith = 21
        integer, parameter :: iufish = 22
        integer, parameter :: iufraw = 23
        integer, parameter :: iufred = 24
      end module ius


      module dims
        integer, parameter :: mx = 256
        integer, parameter :: my = 256
        integer, parameter :: mc = 4
        integer, parameter :: mp = 256
        integer, parameter :: mq = 4
        integer, parameter :: mr = 24
        integer, parameter :: ms = 16384
        integer, parameter :: mt = 128
        integer, parameter :: mu = 4096
        integer, parameter :: mv = 256
        integer, parameter :: mw = 128
        integer, parameter :: mbytes = 2880
        integer, parameter :: mlines = 36
        integer, parameter :: mdata = 720
      end module dims


      module modes
        integer, parameter :: MUNKNOWN = 0
        integer, parameter :: MSTARE = 1
        integer, parameter :: MFLAT = 2
        integer, parameter :: MNOD = 3
        integer, parameter :: MCHOP = 4
        integer, parameter :: MCHOPNOD = 5
        integer, parameter :: MSCAN = 6
        integer, parameter :: MMAP = 7
        integer, parameter :: MHIMED = 1
        integer, parameter :: MHILOW = 2
        integer, parameter :: MMED = 3
        integer, parameter :: MLOW = 4
        integer, parameter :: MCAMERA = 5
        integer, parameter :: MNONE = 1
        integer, parameter :: MBLK = 2
        integer, parameter :: MSKY = 3
        integer, parameter :: MSHINY = 4
        integer, parameter :: MBLKSKY = 5
        integer, parameter :: MBLKSHINY = 6
        integer, parameter :: MBLKOBJ = 7
        integer, parameter :: MOBJ = 8
        integer, parameter :: MCELL = 9
        integer, parameter :: MOLD = 10
        integer, parameter :: MUNWT = 2
        integer, parameter :: MNODWT = 3
        integer, parameter :: MCONWT = 4
        integer, parameter :: MINTWT = 5
        integer, parameter :: MCONINT = 6
        integer, parameter :: MAVGWT = 7
        integer, parameter :: ITNONE = 0
        integer, parameter :: ITHEAD = 1
        integer, parameter :: ITDATA = 2
      end module modes


      module consts
        real, parameter :: DEGRAD = 57.2958
        real, parameter :: PI = 3.14159
        real, parameter :: TPI = 6.28319
        real, parameter :: TINY = 1.e-06
        real, parameter :: HUGE = 1.0e+31
        real, parameter :: CLIGHT = 2.998e+05
        real, parameter :: hc2 = 5.955e-06
        real, parameter :: hck = 1.4388
      end module consts

      module paths
        character(8), parameter :: pgdev = '/xwindow'
        character(32), parameter :: atmodir = 'atmodir not set'
      end module paths
