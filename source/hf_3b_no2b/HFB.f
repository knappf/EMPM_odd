      program Phonon_EDF

       USE technical
       USE math

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

        call read_input

        call dimm(id,noscmax)  !calculates the dimension 'id' for no. of shells 'noscmax'

        call initialize_basis
        call occupations

        call kinetic
        write(*,*)'Loading interaction ...'
c        call interaction
        call interaction_binary
        write(*,*)'Interaction loaded'

        call transition 

        write(*,*)'HFB iteration starts'
        call hfb_iteration
        write(*,*)'HFB iteration ends'

        call transf_interaction

        call TDA

c        call pnTDApp

c        call pnTDAhh

c        call pnTDAph

c        call pnTDAhp

c        call ppTDApp

c        call nnTDApp

c        call TDAtz

c        call deallocate_all

      end

