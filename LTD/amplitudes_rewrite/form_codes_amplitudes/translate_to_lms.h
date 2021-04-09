#procedure introduce-lms
* REPLACE EVERYTHING APART FROM LOOP-MOMENTA BY lm's
* the argument environments nesting is needed for replacing the props
* replace scalar products involving loop-momenta    
    #do i=`NFINALMOMENTA'+1,`NCUTMOMENTA'
        id c`i'.p1? = penergy(c`i')*penergy(p1) - spatial(c`i',p1);
        symmetrize spatial;
        argument;
            id c`i'.p1? = penergy(c`i')*penergy(p1) - spatial(c`i',p1);
            symmetrize spatial;
        endargument;
    #enddo
    .sort:expand-sps-loop-mom;
***************** TRANSLATION OF FUNCTIONS INVOLVING ONLY MOMENTA (as in standard LTD)
* lm-replacements: we do not replace loop-energies, since they are needed for the PF-routine
* Convert the dot products and energies to a symbol
    #$MAXK = `NCUTMOMENTA';
    #$MAXP = `NINITIALMOMENTA';
    #$OFFSET = 0;
    #do i=1,`$MAXP'
        id penergy(p`i') = lm`$OFFSET';
        argument;
            id penergy(p`i') = lm`$OFFSET';
            argument;
                id penergy(p`i') = lm`$OFFSET';
            endargument;
        endargument;

        #$OFFSET = $OFFSET + 1;
        #do j=`i',`$MAXP'
            id p`i'.p`j' = lm`$OFFSET';
            argument;
                id p`i'.p`j' = lm`$OFFSET';
                argument;
                    id p`i'.p`j' = lm`$OFFSET';
                endargument;
            endargument;
            #$OFFSET = $OFFSET + 1;
            id spatial(p`i', p`j') = lm`$OFFSET';
            argument;
                id spatial(p`i', p`j') = lm`$OFFSET';
                argument;
                    id spatial(p`i', p`j') = lm`$OFFSET';
                endargument;
            endargument;
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo

    #do i=1,`$MAXK'
        #if (`i'<= `NFINALMOMENTA');
            id penergy(c`i') = lm`$OFFSET';
            argument;
                id penergy(c`i') = lm`$OFFSET';
                argument;
                    id penergy(c`i') = lm`$OFFSET';
                endargument;
            endargument;
        #endif
        #$OFFSET = $OFFSET + 1;
        #do j=1,`$MAXP'
* there sould not exist any, which involve loop-momenta, because sp's involving loop-momenta are expanded beforehand
            id c`i'.p`j' = lm`$OFFSET';
            argument;
                id c`i'.p`j' = lm`$OFFSET';
                argument;
                    id c`i'.p`j' = lm`$OFFSET';
                endargument;            
            endargument;
            #$OFFSET = $OFFSET + 1;
            id spatial(p`j', c`i') = lm`$OFFSET';
            argument;
                id spatial(p`j', c`i') = lm`$OFFSET';
                argument;
                    id spatial(p`j', c`i') = lm`$OFFSET';
                endargument;
            endargument;
            #$OFFSET = $OFFSET + 1;
        #enddo

        #do j=`i',`$MAXK'
            id c`i'.c`j' = lm`$OFFSET';
            argument;
                id c`i'.c`j' = lm`$OFFSET';
                argument;
                    id c`i'.c`j' = lm`$OFFSET';
                endargument;
            endargument;
            #$OFFSET = $OFFSET + 1;
            id spatial(c`i', c`j') = lm`$OFFSET';
            argument;
                id spatial(c`i', c`j') = lm`$OFFSET';
                argument;
                    id spatial(c`i', c`j') = lm`$OFFSET';
                endargument;
            endargument;
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo
********************** TRANSLATION OF ALL ADDITONAL OBJECTS ONLY EXISTING IN AMPLITUDES (spinors, spatial components etc)
* for polarized cross-sections
    #$MAXEPS =  `NPOL';
    #$MAXCEPS = `NCPOL';
    #$MAXV =`NSPINV';
    #$MAXVBAR=  `NSPINVBAR';
    #$MAXU  = `NSPINU';
    #$MAXUBAR = `NSPINUBAR';
* spatial components of momenta    
    #do i=1,`$MAXP'
        #do j =1,3
            id spatialComp(p`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo       
    #enddo
    #do i=1,`$MAXK'
        #do j =1,3
            id spatialComp(c`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo
* polarization vectors
    #do i=1,`$MAXEPS'
        id penergy(eps`i') = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
        #do j =1,3
            id spatialComp(eps`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
        #do j=`i', `$MAXEPS'
            id eps`i'.eps`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
        #enddo
        #do j=1, `$MAXK'
            id eps`i'.c`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
            id spatial(c`j', eps`i') = lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo       
        #do j=1, `$MAXP'
            id eps`i'.p`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
        #enddo        
        #do j=1, `$MAXCEPS'
            id eps`i'.ceps`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
        #enddo
    #enddo

    #do i=1,`$MAXCEPS'
        id penergy(ceps`i') =   lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
        #do j =1,3
            id spatialComp(ceps`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
        #do j=`i', `$MAXCEPS'
            id ceps`i'.ceps`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
        #enddo
        #do j=1, `$MAXK'
            id ceps`i'.c`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
            id spatial(c`j', ceps`j') = lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo       
        #do j=1, `$MAXP'
            id ceps`i'.p`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
        #enddo    
    #enddo
* spinors    
    #do i=1,`$MAXV'
        id penergy(sV`i') =   lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
        #do j =1,3
            id spatialComp(sV`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo
    #do i=1,`$MAXVBAR'
        id penergy(sVbar`i') = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
        #do j =1,3
            id spatialComp(sVbar`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo
    #do i=1,`$MAXU'
        id penergy(sU`i')  = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
        #do j =1,3
            id spatialComp(sU`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo
    #do i=1,`$MAXUBAR'
        id penergy(sUbar`i') = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
        #do j =1,3
            id spatialComp(sUbar`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo
    argument prop;
        id sqrt(x?) = (x)^(1/2);
    endargument;
#endprocedure