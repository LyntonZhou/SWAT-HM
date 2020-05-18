    real function hmlresuspension(d)

    !!     ~ ~ ~ PURPOSE ~ ~ ~
    !!    this subroutine computes the resuspension velpcity
    !
    !!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
    !!    name          |units         |definition
    !!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !!    ch_l2(:)      |km            |length of main channel
    !!    ch_w(2,:)     |m             |average width of main channel
    !!    ch_s(2,:)     |m/m           |average slope of main channel
    !!    phi(6,:)      |m             |bottom width of main channel
    !!    chside(:)     |none          |change in horizontal distance per unit
    !!                                 |change in vertical distance on channel side
    !!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    !!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
    !!    name        |units         |definition
    !!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !!    d           |m             | particle dimeter
    !!    vr          |m/s           | resuspension vlocity
    !!    rh          |m             | hydraulic radius
    !!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    use parm


    integer :: jrch
    real :: c1, c2, c3, c4
    real :: c, p
    real :: ps
    real :: d, n, slope, depth, rh, vr, vc
    real :: shearstress, shearstress_cn, shearstress_c
    real :: a !! adjust factor
    real :: zlf

    !!  initialize
    c1 = 8.5*10**-16.
    c2 = 9.07
    c3 = 8.64*10**3
    c4 = 7.
    n = 1.
    ps = 1.26      !! particle density
    !d = 10*10**-6. !! particle dimeter 10 um
    a = 25.
    jrch = 0
    jrch = inum1
    vc = vel_chan(jrch)

    if (vc > 5.) vc = 5.

    depth = 0.
    if (rchdep < 0.01) then
        depth = .01
    else
        depth = rchdep
    endif
    ! added by Zhou for debug
    if (depth > 2.5) then
        zlf = 1.
    endif
    if (vc > 4.99) then
        zlf = 1.
    endif

    ! depth = 1.0
    c = chside(jrch)
    slope = ch_s(2,irch)
    rcharea = (phi(6,jrch) + c * depth) * depth
    p = phi(6,jrch) + 2. * depth * Sqrt(1. + c * c)
    rh = rcharea / p
    !shearstress = 1000*9.8*rh*slope/a
    shearstress = 1000*0.003*vc**2
    !! noncohesive sediment critial sheer stress
    shearstress_cn = 414*d
    !! cohesive sediment critial sheer stress
    shearstress_c = shearstress_cn*(1+c1*exp(c2*ps)/(d*10**6.)**2 + c3/(c4*(d*10**6.)))

    if (shearstress < shearstress_cn) then
        vr = 0.
    else
        vr = ((shearstress-shearstress_cn)/(shearstress_c-shearstress_cn))**n * 10**-6.
    endif

    vr = vr*3600*24 !! m/s to m/day

    !if (vr > 10.) then
    !    vr = 10.
    !end if

    hmlresuspension = vr

    return
    end
