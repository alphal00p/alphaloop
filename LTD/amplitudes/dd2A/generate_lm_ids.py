generate_lm_ids(n_cut_momenta, n_initial_momenta, n_final_momenta):
    offset = 0
    id_list = []
    for i in range(n_initial_momenta)
        id_list += ["penergy(p%d) = lm%d" % (i+1, offset)]
        offset += 1
        for j in range(i, n_initial_momenta)
            id_list["p%d.p%d = lm%d"%(i+1,j+1,offset)]
            offset += 1
            id_list["spatial(p%d,p%d) = lm%d"%(i+1,j+1,offset)]
            offset += 1

    for i in range(n_cut_momenta)
        if i+1 <= n_final_momenta: 
            id_list["penergy(c%d) = lm%d"%(i+1,offset)]
        offset += 1
        for j in range(n_initial_momenta)
            #there sould not exist any, which involve loop-momenta, 
            # because sp's involving loop-momenta are expanded beforehand
            id_list["c%d.p%d = lm%d"%(i+1,j+1,offset)]
            offset += 1
            id_list["spatial(c%d,p%d) = lm%d"%(i+1,j+1,offset)]
            offset += 1

        for j in range(i, n_cut_momenta)
            id_list["c%d.c%d = lm%d"%(i+1,j+1,offset)]
            offset += 1
            id_list["spatial(c%d,c%d) = lm%d"%(i+1,j+1,offset)]
            offset += 1
    #translation of all additonal objects only existing in amplitudes(spinors, spatial components etc)
* for polarized cross-sections
    # $MAXEPS =  `NPOL';
    # $MAXCEPS = `NCPOL';
    # $MAXV =`NSPINV';
    # $MAXVBAR=  `NSPINVBAR';
    # $MAXU  = `NSPINU';
    # $MAXUBAR = `NSPINUBAR';
* spatial components of momenta
    # do i=1,`n_initial_momenta'
        # do j =1,3
            id spatialComp(p`i',`j')=lm`$OFFSET';
            offset += 1
        # enddo
    # enddo
    # do i=1,`n_cut_momenta'
        # do j =1,3
            id spatialComp(c`i',`j')=lm`$OFFSET';
            offset += 1
        # enddo
    # enddo
* polarization vectors
    # do i=1,`$MAXEPS'
        id penergy(eps`i') = lm`$OFFSET';
        offset += 1

        # do j =1,3
            id spatialComp(eps`i',`j')=lm`$OFFSET';
            offset += 1
        # enddo

        # do j=`i', `$MAXEPS'
            id eps`i'.eps`j'=lm`$OFFSET';
            # $OFFSET = $OFFSET+1;
        # enddo

        # do j=1, `n_cut_momenta'
            id eps`i'.c`j'=lm`$OFFSET';
            # $OFFSET = $OFFSET+1;

            id spatial(c`j', eps`i')=lm`$OFFSET';
            offset += 1
        # enddo

        # do j=1, `n_initial_momenta'
            id eps`i'.p`j'=lm`$OFFSET';
            # $OFFSET = $OFFSET+1;
        # enddo

        # do j=1, `$MAXCEPS'
            id eps`i'.ceps`j'=lm`$OFFSET';
            # $OFFSET = $OFFSET+1;
        # enddo
    # enddo

    # do i=1,`$MAXCEPS'
        id penergy(ceps`i') =   lm`$OFFSET';
        offset += 1

        # do j =1,3
            id spatialComp(ceps`i',`j')=lm`$OFFSET';
            offset += 1
        # enddo
        # do j=`i', `$MAXCEPS'
            id ceps`i'.ceps`j'=lm`$OFFSET';
            # $OFFSET = $OFFSET+1;
        # enddo
        # do j=1, `n_cut_momenta'
            id ceps`i'.c`j'=lm`$OFFSET';
            # $OFFSET = $OFFSET+1;

            id spatial(c`j', ceps`i')=lm`$OFFSET';
            offset += 1
        # enddo
        # do j=1, `n_initial_momenta'
            id ceps`i'.p`j'=lm`$OFFSET';
            # $OFFSET = $OFFSET+1;
        # enddo
    # enddo
* spinors
    # do i=1,`$MAXV'
        id penergy(sV`i') =   lm`$OFFSET';

        offset += 1
        # do j =1,3
            id spatialComp(sV`i',`j')=lm`$OFFSET';
            offset += 1
        # enddo
    # enddo

    # do i=1,`$MAXVBAR'
        id penergy(sVbar`i') = lm`$OFFSET';
        offset += 1
        # do j =1,3
            id spatialComp(sVbar`i',`j')=lm`$OFFSET';
            offset += 1
        # enddo
    # enddo
    # do i=1,`$MAXU'
        id penergy(sU`i')  = lm`$OFFSET';
        offset += 1
        # do j =1,3
            id spatialComp(sU`i',`j')=lm`$OFFSET';
            offset += 1
        # enddo
    # enddo
    # do i=1,`$MAXUBAR'
        id penergy(sUbar`i') = lm`$OFFSET';
        offset += 1
        # do j =1,3
            id spatialComp(sUbar`i',`j')=lm`$OFFSET';
            offset += 1
        # enddo
    # enddo

# endprocedure
