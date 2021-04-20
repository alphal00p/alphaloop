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
    Format Mathematica;
    #if (`DEBUGLVL'>0)
        #write<debug_diag_`SGID'.m> "lmRepl`SGID'={"
    #endif
.sort-debug-output;
***************** TRANSLATION OF FUNCTIONS INVOLVING ONLY MOMENTA (as in standard LTD)
* lm-replacements: we do not replace loop-energies, since they are needed for the PF-routine
* Convert the dot products and energies to a symbol
    #$MAXK = `NCUTMOMENTA';
    #$MAXP = `NINITIALMOMENTA';
    #$OFFSET = 0;
    #do i=1,`$MAXP'
        #if (`DEBUGLVL'>0)
            #write<debug_diag_`SGID'.m> "penergy(p`i') -> lm`$OFFSET',"
        #endif
        id penergy(p`i') = lm`$OFFSET';
        argument;
            id penergy(p`i') = lm`$OFFSET';
            argument;
                id penergy(p`i') = lm`$OFFSET';
            endargument;
        endargument;

        #$OFFSET = $OFFSET + 1;
        #do j=`i',`$MAXP'
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "sp(p`i',p`j') -> lm`$OFFSET',"
            #endif
            id p`i'.p`j' = lm`$OFFSET';
            argument;
                id p`i'.p`j' = lm`$OFFSET';
                argument;
                    id p`i'.p`j' = lm`$OFFSET';
                endargument;
            endargument;
            #$OFFSET = $OFFSET + 1;

            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatial(p`i',p`j') -> lm`$OFFSET',"
            #endif
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
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "penergy(c`i') -> lm`$OFFSET',"
            #endif
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
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "sp(c`i',p`j') -> lm`$OFFSET',"
            #endif
            id c`i'.p`j' = lm`$OFFSET';
            argument;
                id c`i'.p`j' = lm`$OFFSET';
                argument;
                    id c`i'.p`j' = lm`$OFFSET';
                endargument;            
            endargument;
            #$OFFSET = $OFFSET + 1;

            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatial(p`j', c`i') -> lm`$OFFSET',"
            #endif
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

            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "sp(c`i',c`j') -> lm`$OFFSET',"
            #endif
            id c`i'.c`j' = lm`$OFFSET';
            argument;
                id c`i'.c`j' = lm`$OFFSET';
                argument;
                    id c`i'.c`j' = lm`$OFFSET';
                endargument;
            endargument;
            #$OFFSET = $OFFSET + 1;

            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatial(c`i',c`j') -> lm`$OFFSET',"
            #endif
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
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatialComp(p`i',`j') -> lm`$OFFSET',"
            #endif
            id spatialComp(p`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo       
    #enddo
    #do i=1,`$MAXK'
        #do j =1,3
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatialComp(c`i',`j') -> lm`$OFFSET',"
            #endif
            id spatialComp(c`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo
* polarization vectors
    #do i=1,`$MAXEPS'
        #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "penergy(eps`i') -> lm`$OFFSET',"
        #endif
        id penergy(eps`i') = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;

        #do j =1,3        
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatialComp(eps`i',`j') -> lm`$OFFSET',"
            #endif
            id spatialComp(eps`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo

        #do j=`i', `$MAXEPS'
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "sp(eps`i',eps`j') -> lm`$OFFSET',"
            #endif
            id eps`i'.eps`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
        #enddo

        #do j=1, `$MAXK'
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "sp(eps`i',c`j') -> lm`$OFFSET',"
            #endif
            id eps`i'.c`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;

            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatial(eps`i',c`j') -> lm`$OFFSET',"
            #endif
            id spatial(c`j', eps`i') = lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo       

        #do j=1, `$MAXP'
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "sp(eps`i',p`j') -> lm`$OFFSET',"
            #endif
            id eps`i'.p`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
        #enddo        

        #do j=1, `$MAXCEPS'
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "sp(eps`i',ceps`j') -> lm`$OFFSET',"
            #endif
            id eps`i'.ceps`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
        #enddo
    #enddo

    #do i=1,`$MAXCEPS'
        #if (`DEBUGLVL'>0)
            #write<debug_diag_`SGID'.m> "penergy(ceps`i') -> lm`$OFFSET',"
        #endif
        id penergy(ceps`i') =   lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;

        #do j =1,3
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatialComp(ceps`i',`j') -> lm`$OFFSET',"
            #endif
            id spatialComp(ceps`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
        #do j=`i', `$MAXCEPS'
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatialComp(ceps`i',ceps`j) -> lm`$OFFSET',"
            #endif
            id ceps`i'.ceps`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
        #enddo
        #do j=1, `$MAXK'
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "sp(ceps`i',c`j') -> lm`$OFFSET',"
            #endif
            id ceps`i'.c`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;

            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatial(c`j', ceps`j') -> lm`$OFFSET',"
            #endif
            id spatial(c`j', ceps`j') = lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo       
        #do j=1, `$MAXP'
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "sp(c`j', p`j') -> lm`$OFFSET',"
            #endif
            id ceps`i'.p`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
        #enddo    
    #enddo
* spinors    
    #do i=1,`$MAXV'
        id penergy(sV`i') =   lm`$OFFSET';
        #if (`DEBUGLVL'>0)
            #write<debug_diag_`SGID'.m> "penergy(sV`i') -> lm`$OFFSET',"
        #endif

        #$OFFSET = $OFFSET + 1;
        #do j =1,3
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatialComp(sV`i',`j') -> lm`$OFFSET',"
            #endif
            id spatialComp(sV`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo

    #do i=1,`$MAXVBAR'
        #if (`DEBUGLVL'>0)
            #write<debug_diag_`SGID'.m> "penergy(sVbar`i') -> lm`$OFFSET',"
        #endif
        id penergy(sVbar`i') = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
        #do j =1,3
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatialComp(sVbar`i',`j') -> lm`$OFFSET',"
            #endif
            id spatialComp(sVbar`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo
    #do i=1,`$MAXU'
        #if (`DEBUGLVL'>0)
            #write<debug_diag_`SGID'.m> "penergy(sU`i') -> lm`$OFFSET',"
        #endif
        id penergy(sU`i')  = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
        #do j =1,3
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatialComp(sU`i',`j') -> lm`$OFFSET',"
            #endif
            id spatialComp(sU`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo
    #do i=1,`$MAXUBAR'
        #if (`DEBUGLVL'>0)
            #write<debug_diag_`SGID'.m> "penergy(sUbar`i') -> lm`$OFFSET',"
        #endif
        id penergy(sUbar`i') = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
        #do j =1,3
            #if (`DEBUGLVL'>0)
                #write<debug_diag_`SGID'.m> "spatialComp(sUbar`i',`j') -> lm`$OFFSET',"
            #endif
            id spatialComp(sUbar`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo

    #if (`DEBUGLVL'>0)
        #write<debug_diag_`SGID'.m> "1->1}"
    #endif
#endprocedure