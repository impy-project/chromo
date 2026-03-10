//
//  This file is part of EPOS4 modified for EPOS-LHCR
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include "ak.h"

using namespace std ;


//###################################################################################################
//################################ maxsize object -> dimensions #####################################
//###################################################################################################

Mudiar<int> *maxsize ;
extern "C" void maxsize_create_(int* n1)      { maxsize = new Mudiar<int>(*n1) ;  }
extern "C" void maxsize_destroy_(void)        { delete maxsize ; }
extern "C" void maxsize_set_(int* i1,int* val){ maxsize->set(*i1-1, *val) ;  }
extern "C" void maxsize_get_(int* i1,int* val){ maxsize->get(*i1-1, *val) ;  }

//###################################################################################################
//################################### event variables object ########################################
//###################################################################################################

//-------------------------------------------------------------------------------------------

Mudiar<float> *eventvari ;
extern "C" void eventvaricreate_(int* n1)         { eventvari = new Mudiar<float>(*n1) ;  }
extern "C" void eventvaridestroy_(void)           { delete eventvari ; }
extern "C" void eventvariset_(int* i1,float* val) { eventvari->set(*i1-1, *val) ;  }
extern "C" void eventvariget_(int* i1,float* val) { eventvari->get(*i1-1, *val) ;  }

//###################################################################################################
//###################################################################################################

//-------------------------------------------------------------------------------------------

Mudiar<int> *lspecs ;
extern "C" void lspecscreate_(int* n1,int* n2)      { lspecs = new Mudiar<int>(*n1,*n2) ;  }
extern "C" void lspecsdestroy_(void)                { delete lspecs ; }
extern "C" void lspecsset_(int* i1,int* i2,int* val){ lspecs->set(*i1-1,*i2-1, *val) ;  }
extern "C" void lspecsget_(int* i1,int* i2,int* val){ lspecs->get(*i1-1,*i2-1, *val) ;  }
extern "C" void lspecsincrement_(int* i1,int* i2)   { lspecs->increment(*i1-1,*i2-1) ;  }

Mudiar<float> *wgtpairst ;
extern "C" void wgtpairstcreate_(int* n1)      { wgtpairst = new Mudiar<float>(*n1) ;  }
extern "C" void wgtpairstdestroy_(void)                { delete wgtpairst ; }
extern "C" void wgtpairstset_(int* i1,float* val){ wgtpairst->set(*i1-1, *val) ;  }
extern "C" void wgtpairstget_(int* i1,float* val){ wgtpairst->get(*i1-1, *val) ;  }

Mudiar<int> *idpairst ;
extern "C" void idpairstcreate_(int* n1,int* n2)      { idpairst = new Mudiar<int>(*n1,*n2) ;  }
extern "C" void idpairstdestroy_(void)                        { delete idpairst ; }
extern "C" void idpairstset_(int* i1,int* i2,int* val){ idpairst->set(*i1-1,*i2-1, *val) ;  }
extern "C" void idpairstget_(int* i1,int* i2,int* val){ idpairst->get(*i1-1,*i2-1, *val) ;  }

Mudiar<int> *lkfok ;
extern "C" void lkfokcreate_(int* n1, int* n2, int* n3, int* n4, int* n5)       { lkfok = new Mudiar<int>(*n1, *n2, *n3, *n4, *n5) ; }
extern "C" void lkfokdestroy_(void)                                             { delete lkfok ; }
extern "C" void lkfokset_(int* i1, int* i2, int* i3, int* i4, int* i5, int* val){ lkfok->set(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void lkfokget_(int* i1, int* i2, int* i3, int* i4, int* i5, int* val){ lkfok->get(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void lkfokincrement_(int* i1, int* i2, int* i3, int* i4, int* i5)    { lkfok->increment(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1) ; }

//-------------------------------------------------------------------------------------------

Mudiar<float> *emuc ;
extern "C" void createemuc_(int* n1, int* n2, int* n3, int* n4, int* n5)          { emuc = new Mudiar<float>(*n1, *n2, *n3, *n4, *n5) ; }
extern "C" void destroyemuc_(void)                                                { delete emuc ; }
extern "C" void emucset_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { emuc->set(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void emucget_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { emuc->get(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
 
//-------------------------------------------------------------------------------------------

Mudiar<float> *velio ;
extern "C" void createvelio_(int* n1, int* n2, int* n3, int* n4, int* n5)          { velio = new Mudiar<float>(*n1, *n2, *n3, *n4, *n5) ; }
extern "C" void destroyvelio_(void)                                                { delete velio ; }
extern "C" void velioset_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { velio->set(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void velioget_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { velio->get(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }

Mudiar<float> *bario ;
extern "C" void createbario_(int* n1, int* n2, int* n3, int* n4, int* n5)          { bario = new Mudiar<float>(*n1, *n2, *n3, *n4, *n5) ; }
extern "C" void destroybario_(void)                                                { delete bario ; }
extern "C" void barioset_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { bario->set(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void barioget_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { bario->get(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }

Mudiar<float> *epsio ;
extern "C" void createepsio_(int* n1, int* n2, int* n3, int* n4, int* n5)          { epsio = new Mudiar<float>(*n1, *n2, *n3, *n4, *n5) ; }
extern "C" void destroyepsio_(void)                                                { delete epsio ; }
extern "C" void epsioset_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { epsio->set(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void epsioget_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { epsio->get(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }

Mudiar<float> *emuio ;
extern "C" void createemuio_(int* n1, int* n2, int* n3, int* n4, int* n5)          { emuio = new Mudiar<float>(*n1, *n2, *n3, *n4, *n5) ; }
extern "C" void destroyemuio_(void)                                                { delete emuio ; }
extern "C" void emuioset_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { emuio->set(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }
extern "C" void emuioget_(int* i1, int* i2, int* i3, int* i4, int* i5, float* val) { emuio->get(*i1-1, *i2-1, *i3-1, *i4-1, *i5-1, *val) ; }

