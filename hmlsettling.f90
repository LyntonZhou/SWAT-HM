    real function hmlsettling(d)

    !!     ~ ~ ~ PURPOSE ~ ~ ~
    !!    this subroutine computes the settling velpcity

    !!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
    !!    name        |units         |definition
    !!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !!    d           |m             | particle dimeter
    !!    vs          |m/s           | settling vlocity
    !!    mu          |g/cm/s        ! viscosity
    !!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    use parm

    integer :: jrch
    real :: ps, pw
    real :: d, mu, vs
    real :: zlf !! debug

    !!  initialize
    ps = 2.65          !! g/cm3
    pw = 1.0           !! g/cm3
    mu = 1.002/100     !! g/cm/s 
    !d = 10/10000       !! cm
    
    d = d *100  ! m to cm
    
    vs = 980*d**2*(ps-pw)/(18*mu) ! cm/s
    vs = vs/100*86400            ! m/day
    
    if (vs > 5.) then
        zlf = 1.
    end if

    hmlsettling = vs
    return
    end
