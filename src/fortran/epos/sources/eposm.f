      subroutine NumberModel(cmodel,model)
      character cmodel*21,cmodel2*21
      cmodel2=cmodel
      model2=model
      stop'   ***** This program can only run EPOS *****      '
      end

      subroutine IniModel(model)
      model2=model
      end

      subroutine IniEvtModel
      end

      subroutine emsaaaModel(model,id,iret)
      model2=model
      id2=id
      iret2=iret
      end

      function crseModel(model,ekin,maproj,matarg,idtarg)
      crseModel=0.
      model2=model
      ekin2=ekin
      maproj2=maproj
      matarg2=matarg
      idtarg2=idtarg
      end

      subroutine crseaaModel(sigt,sigi,sigc,sigel)
      sigt=0.
      sigi=0.
      sigc=0.
      sigel=0.
      end

      subroutine m2XXFZ( a,b)
      double precision a, b(2),c,d(2)
      c=a
      d(1)=b(1)
      d(2)=b(2)
      
      end

      subroutine m3SIGMA(ek,idpro,idtar,latar,matar,sigi,sige)
      sige=0.
      sigi=0.
      ek2=ek
      idpro2=idpro
      matar2=matar
      idtar2=idtar
      latar2=latar
      end

      subroutine m6SIGMA(icl,engy,stot,sela,sine,sdifr,slela,Rho)
      icl2=icl
      engy2=engy
      stot=0.
      sela=0.
      sine=0.
      sdifr=0.
      slela=0.
      Rho=0.
      end

      subroutine m7SIGMA(stot,scut,sine,slela)
      stot=0.
      scut=0.
      sine=0.
      slela=0.
      end

      subroutine m8SIGMA(stot,scut,sine,sela,slela,ssd)
      stot=0.
      scut=0.
      sine=0.
      sela=0.
      slela=0.
      ssd=0.
      end

      subroutine m9SIGMA(stot,sine,sela)
      stot=0.
      sine=0.
      sela=0.
      end

      subroutine decaymod(ip)
      idum=ip
      end
