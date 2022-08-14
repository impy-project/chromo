c $Id: options.f,v 1.9 2001/04/06 21:48:16 weber Exp $
c... law: include file (only) for global parameters & options
      integer numcto,numctp,maxstables
      parameter(numcto=400)
      parameter(numctp=400)
      parameter(maxstables=20)
c...
      integer   CTOption(numcto)
      character ctodc(numcto)*2
c...
      real*8    CTParam(numctp)
      character ctpdc(numctp)*2

      integer nstable
      integer stabvec(maxstables)

      logical bf13,bf14,bf15,bf16,bf17,bf18,bf19,bf20,fixedseed
      common /options/CTOption,CTParam
      common /optstrings/ctodc,ctpdc
      common /loptions/fixedseed,bf13,bf14,bf15,bf16,bf17,bf18,
     .     bf19,bf20
      common /stables/nstable,stabvec
