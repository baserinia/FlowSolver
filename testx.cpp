// CMeshAdap: meshadap.cpp
// Mesh object implementation
//
// Author: A.R. Baserinia
// Date: September 17, 2004
// Update: September 20, 2004
// Note: This module is only applicable to triangular grids, 
//-----------------------------------------------------------------------------

#include <vector>
#include <algorithm>
#include "meshadap.h"
#include "intface.h"
#include "vector2D.h"

#define PI	3.1415926535

using namespace std;
using namespace NBCFD;


//-------------------------------------
// Split Face						  
//-------------------------------------
void CMeshAdap::SplitFace( CFace* pFace )
{
	int i;
	int j;

	if (pFace->GetCellNum() == 0) return;

	CVertex* lpVertex = new CVertex( pFace->GetFMP() );
	lpVertex->SetType( pFace->GetCellNum() == 2? 0 : 1 );  // 0 for interior and 1 for boundary
	m_Verts.push_back( lpVertex );
	
	if ( pFace->GetCellNum() == 2 ) {	// internal face

		// create temporary cells
		int liNC =	static_cast<CCell*>( pFace->GetCell( 0 ) )->GetFaceNum() +
					static_cast<CCell*>( pFace->GetCell( 1 ) )->GetFaceNum() - 2;
		vector<CCell*> lTmpCell( liNC );
		lTmpCell[0] = static_cast<CCell*>( pFace->GetCell( 0 ) );
		lTmpCell[1] = static_cast<CCell*>( pFace->GetCell( 1 ) );
		for (i=2; i<liNC; i++) {
			lTmpCell[i] = new CCell( );
			m_Cells.push_back( lTmpCell[i] );
		}

		// create temporary faces
		int liNF =	static_cast<CCell*>( pFace->GetCell( 0 ) )->GetVertNum() +
					static_cast<CCell*>( pFace->GetCell( 1 ) )->GetVertNum() - 2;
		vector<CFace*> lTmpFace( liNF );
		lTmpFace[0] = pFace;
		for (i=1; i<liNF; i++) {
			lTmpFace[i] = new CIntFace();
			m_Faces.push_back( lTmpFace[i] );
		}

		
		vector<CFace*> lTmpBndF( liNF );	// cavity boundary faces
		j = 0;
		for (i=0; i<lTmpCell[0]->GetFaceNum(); i++) {
			if ( lTmpCell[0]->GetFace( i ) != pFace ) {
				lTmpBndF[j] =  static_cast<CFace*>( lTmpCell[0]->GetFace( i ) ) ;
				j++;
			}
		}
		for (i=0; i<lTmpCell[1]->GetFaceNum(); i++) {
			if ( lTmpCell[1]->GetFace( i ) != pFace ) {
				lTmpBndF[j] =  static_cast<CFace*>( lTmpCell[1]->GetFace( i ) ) ;
				j++;
			}
		}


		vector<CVertex*> lTmpVert( liNF );	// cavity boundary vertices
		j = 0;
		for (i=0; i<liNF; i++) {
			if ( lTmpVert.end() ==	find(	lTmpVert.begin(), 
											lTmpVert.end(), 
											lTmpBndF[i]->GetVert(0) ) ) {
				 lTmpVert[j] = static_cast<CVertex*>( lTmpBndF[i]->GetVert(0) );
				 j++;
			}
			if ( lTmpVert.end() ==	find(	lTmpVert.begin(), 
											lTmpVert.end(), 
											lTmpBndF[i]->GetVert(1) ) ) {
				 lTmpVert[j] = static_cast<CVertex*>( lTmpBndF[i]->GetVert(1) ) ;
				 j++;
			}
		}

		CVertex* pV0 = lTmpVert[0];
		CVertex* pV1 = lTmpVert[1];
		CVertex* pV2 = lTmpVert[2];
		CVertex* pV3 = lTmpVert[3];

		for (i=0; i<liNF; i++) {
			if ( lTmpVert[i]->IsNbrFace( pFace ) ) {
				if ( lTmpVert[i] == pFace->GetVert(0) )
					lTmpVert[i]->RemoveVert( pFace->GetVert(1) );
				if ( lTmpVert[i] == pFace->GetVert(1) )
					lTmpVert[i]->RemoveVert( pFace->GetVert(0) );
				lTmpVert[i]->RemoveFace( pFace );
				lTmpVert[i]->RemoveCell( pFace->GetCell(0) );
				lTmpVert[i]->RemoveCell( pFace->GetCell(1) );
			} else {
				if ( lTmpVert[i]->IsNbrCell( pFace->GetCell(0) ) )
					lTmpVert[i]->RemoveCell( pFace->GetCell(0) );
				if ( lTmpVert[i]->IsNbrCell( pFace->GetCell(1) ) )
					lTmpVert[i]->RemoveCell( pFace->GetCell(1) );
			}
		}
		pFace->ClearNbrAll();

		for (i=0; i<liNF; i++) {
			if ( ( lTmpBndF[i]->GetCell(0) == lTmpCell[0] ) ||
				( lTmpBndF[i]->GetCell(0) == lTmpCell[1] ) ) {
				if ( lTmpBndF[i]->GetCellNum() > 1) 
					CAST(CCell*)( lTmpBndF[i]->GetCell(1) )->RemoveCell( lTmpBndF[i]->GetCell(0) );
				lTmpBndF[i]->SetCell( 0, lTmpCell[i] );
			}
			if ( ( lTmpBndF[i]->GetCellNum() > 1) &&
				 ( ( lTmpBndF[i]->GetCell(1) == lTmpCell[0] ) ||
				 ( lTmpBndF[i]->GetCell(1) == lTmpCell[1] ) ) ) {
				CAST(CCell*)( lTmpBndF[i]->GetCell(0) )->RemoveCell( lTmpBndF[i]->GetCell(1) );
				lTmpBndF[i]->SetCell( 1, lTmpCell[i] );
			}
		}

		for (i=0; i<liNF; i++) {
			lTmpCell[i]->ClearNbrAll();
			
			lTmpCell[i]->AddVert( static_cast<CVertex*>( lTmpBndF[i]->GetVert(0) ) );
			lTmpCell[i]->AddVert( static_cast<CVertex*>( lTmpBndF[i]->GetVert(1) ) );
			lTmpCell[i]->AddVert( lpVertex );

			lTmpCell[i]->AddFace( lTmpBndF[i] );
		}

		for (i=0; i<liNF; i++) {
			lTmpFace[i]->ClearNbrVerts();
			lTmpFace[i]->AddVert( lpVertex );
			lTmpFace[i]->AddVert( lTmpVert[i] );

			for (j=0; j<liNF; j++) 
				if ( lTmpVert[i]->IsNbrFace( lTmpBndF[j] ) ) {
					lTmpFace[i]->AddCell( lTmpCell[j] );
					lTmpCell[j]->AddFace( lTmpFace[i] );
				}
		}

		for (i=0; i<liNF; i++) {
			static_cast<CCell*>( lTmpFace[i]->GetCell(0) )->AddCell( lTmpFace[i]->GetCell(1) );
			static_cast<CCell*>( lTmpFace[i]->GetCell(1) )->AddCell( lTmpFace[i]->GetCell(0) );

			if (  lTmpBndF[i]->GetCellNum() > 1 ) { 
				CAST(CCell*)( lTmpBndF[i]->GetCell(0) )->AddCell( lTmpBndF[i]->GetCell(1) );
				CAST(CCell*)( lTmpBndF[i]->GetCell(1) )->AddCell( lTmpBndF[i]->GetCell(0) );
			}
		}

		for (i=0; i<liNF; i++) {
			lTmpVert[i]->AddVert( lpVertex );
			lpVertex->AddVert( lTmpVert[i] );
			lTmpVert[i]->AddFace( lTmpFace[i] );
			lpVertex->AddFace( lTmpFace[i] );
			lTmpVert[i]->AddCell( lTmpFace[i]->GetCell(0) );
			lTmpVert[i]->AddCell( lTmpFace[i]->GetCell(1) );
			lpVertex->AddCell( lTmpCell[i] );
		}

		for ( i=0; i<liNC; i++)
			lTmpCell[i]->UpdateGeom();

		for ( i=0; i<liNC; i++)
			lTmpBndF[i]->UpdateGeom();

		for ( i=0; i<liNF; i++)
			lTmpFace[i]->UpdateGeom();
	} 
	
	else {  // BOUNDARY FACE

		// create temporary cells
		int liNC =	static_cast<CCell*>( pFace->GetCell( 0 ) )->GetFaceNum()-1;
		vector<CCell*> lTmpCell( liNC );
		lTmpCell[0] = static_cast<CCell*>( pFace->GetCell( 0 ) );
		for (i=1; i<liNC; i++) {
			lTmpCell[i] = new CCell( );
			m_Cells.push_back( lTmpCell[i] );
		}

		// create temporary faces
		int liNF =	static_cast<CCell*>( pFace->GetCell( 0 ) )->GetVertNum();
		vector<CFace*> lTmpFace( liNF );
		lTmpFace[0] = pFace;
		for (i=1; i<liNF; i++) {
			lTmpFace[i] = new CIntFace();
			m_Faces.push_back( lTmpFace[i] );
		}

		
		vector<CFace*> lTmpBndF( liNC );	// cavity boundary faces
		j = 0;
		for (i=0; i<lTmpCell[0]->GetFaceNum(); i++) {
			if ( lTmpCell[0]->GetFace( i ) != pFace ) {
				lTmpBndF[j] =  static_cast<CFace*>( lTmpCell[0]->GetFace( i ) ) ;
				j++;
			}
		}

		vector<CVertex*> lTmpVert( liNF );	// cavity boundary vertices
		for (i=0; i<liNF; i++) {
			 lTmpVert[i] = static_cast<CVertex*>( lTmpCell[0]->GetVert(i) );
		}

		for (i=0; i<liNF; i++) {
			if ( lTmpVert[i]->IsNbrFace( pFace ) ) {
				if ( lTmpVert[i] == pFace->GetVert(0) )
					lTmpVert[i]->RemoveVert( pFace->GetVert(1) );
				if ( lTmpVert[i] == pFace->GetVert(1) )
					lTmpVert[i]->RemoveVert( pFace->GetVert(0) );
				lTmpVert[i]->RemoveFace( pFace );
				lTmpVert[i]->RemoveCell( pFace->GetCell(0) );
			} else {
				if ( lTmpVert[i]->IsNbrCell( pFace->GetCell(0) ) )
					lTmpVert[i]->RemoveCell( pFace->GetCell(0) );
			}
		}
		pFace->ClearNbrAll();

		for (i=0; i<liNC; i++) {
			if ( lTmpBndF[i]->GetCell(0) == lTmpCell[0] )  {
				if ( lTmpBndF[i]->GetCellNum() > 1) 
					CAST(CCell*)( lTmpBndF[i]->GetCell(1) )->RemoveCell( lTmpBndF[i]->GetCell(0) );
				lTmpBndF[i]->SetCell( 0, lTmpCell[i] );
			}
			if ( ( lTmpBndF[i]->GetCellNum() > 1) &&
				 ( lTmpBndF[i]->GetCell(1) == lTmpCell[0] ) ) {
				if ( lTmpBndF[i]->GetCellNum() > 1) 
					CAST(CCell*)( lTmpBndF[i]->GetCell(0) )->RemoveCell( lTmpBndF[i]->GetCell(1) );
				lTmpBndF[i]->SetCell( 1, lTmpCell[i] );
			}
		}

		for (i=0; i<liNC; i++) {
			lTmpCell[i]->ClearNbrAll();
			
			lTmpCell[i]->AddVert( static_cast<CVertex*>( lTmpBndF[i]->GetVert(0) ) );
			lTmpCell[i]->AddVert( static_cast<CVertex*>( lTmpBndF[i]->GetVert(1) ) );
			lTmpCell[i]->AddVert( lpVertex );

			lTmpCell[i]->AddFace( lTmpBndF[i] );
		}

		for (i=0; i<liNF; i++) {

			lTmpFace[i]->ClearNbrVerts();
			lTmpFace[i]->AddVert( lpVertex );
			lTmpFace[i]->AddVert( lTmpVert[i] );

			for (j=0; j<liNC; j++) 
				if ( lTmpVert[i]->IsNbrFace( lTmpBndF[j] ) ) {
					lTmpFace[i]->AddCell( lTmpCell[j] );
					lTmpCell[j]->AddFace( lTmpFace[i] );
				}
		}

		for (i=0; i<liNF; i++) 
			if ( lTmpFace[i]->GetCellNum() > 1) {
				CAST(CCell*)( lTmpFace[i]->GetCell(0) )->AddCell( lTmpFace[i]->GetCell(1) );
				CAST(CCell*)( lTmpFace[i]->GetCell(1) )->AddCell( lTmpFace[i]->GetCell(0) );
			}

		for (i=0; i<liNC; i++) 
			if (  lTmpBndF[i]->GetCellNum() > 1 ) { 
				CAST(CCell*)( lTmpBndF[i]->GetCell(0) )->AddCell( lTmpBndF[i]->GetCell(1) );
				CAST(CCell*)( lTmpBndF[i]->GetCell(1) )->AddCell( lTmpBndF[i]->GetCell(0) );
			}
		

		for (i=0; i<liNF; i++) {
			lTmpVert[i]->AddVert( lpVertex );
			lpVertex->AddVert( lTmpVert[i] );
			lTmpVert[i]->AddFace( lTmpFace[i] );
			lpVertex->AddFace( lTmpFace[i] );
			lTmpVert[i]->AddCell( lTmpFace[i]->GetCell(0) );
			if ( lTmpFace[i]->GetCellNum() > 1)
				lTmpVert[i]->AddCell( lTmpFace[i]->GetCell(1) );
		}
		lpVertex->AddCell( lTmpCell[0] );
		lpVertex->AddCell( lTmpCell[1] );

		for ( i=0; i<liNC; i++)
			lTmpCell[i]->UpdateGeom();

		for ( i=0; i<liNC; i++)
			lTmpBndF[i]->UpdateGeom();

		for ( i=0; i<liNF; i++)
			lTmpFace[i]->UpdateGeom();
	}


}


//-------------------------------------
// Collapse Face						  
//-------------------------------------
list<CFace*>::iterator CMeshAdap::CollapseFace( list<CFace*>::iterator   FIter)
{
	// special treatment for trangles
	int i;
	int j;
	CFace* pFace = *FIter;
	CCell* pCell;
	CVertex* pVert0 = CAST(CVertex*)( pFace->GetVert(0) );
	CVertex* pVert1 = CAST(CVertex*)( pFace->GetVert(1) );
	CVertex* pSwap;
	CFace* pFace0;
	CFace* pFace1;
	CCell* pCell0 = NULL;
	CCell* pCell1 = NULL;
	list<CFace*>::iterator   lFIter = ++FIter;
	
	if ( ( pVert0->GetType() == 2 ) && ( pVert1->GetType() == 2 ) ) 
		return lFIter;


	if ( ( pVert0->GetType() > 0 ) && ( pVert1->GetType() > 0 ) &&
		 ( pFace->GetCellNum() == 2 ) )
		return lFIter;



	// Vert0 is retained and Vert1 is removed
	if ( pVert1->GetType() > pVert0->GetType() ) {
		pSwap = pVert1;
		pVert1 = pVert0;
		pVert0 = pSwap;
	}

	for (i=0; i<pVert0->GetCellNum(); i++) {
		if ( CAST(CCell*)( pVert0->GetCell(i) )->IsValid( pVert0, pFace->GetCollapsePos() ) == false )
			return lFIter;
	}

	for (i=0; i<pVert1->GetCellNum(); i++) {
		if ( CAST(CCell*)( pVert1->GetCell(i) )->IsValid( pVert1, pFace->GetCollapsePos() ) == false )
			return lFIter;
	}


	for (i=0; i < pFace->GetCellNum(); i++) {
		// for triangular mesh
		if ( CAST(CCell*)( pFace->GetCell(i) )->GetFaceNum() == 3 ) {
			pCell0 = NULL;
			pCell1 = NULL;
			pFace0 = NULL;
			pFace1 = NULL;

			pCell = CAST(CCell*)( pFace->GetCell(i) );
			for (j=0; j<pCell->GetFaceNum(); j++) {
				if ( pFace != CAST(CFace*)( pCell->GetFace(j) ) ) {
					if ( CAST(CFace*)( pCell->GetFace(j) )->IsNbrVert( pVert0 ) ) {
						pFace0 = CAST(CFace*)( pCell->GetFace(j) );
					} else {
						pFace1 = CAST(CFace*)( pCell->GetFace(j) );
					}
				}
			}

			
			if ( pFace0->GetCellNum() > 1) {
				if ( pFace0->GetCell(0) == pCell ) {
					pCell0 = CAST(CCell*)( pFace0->GetCell(1) );
				} else {
					pCell0 = CAST(CCell*)( pFace0->GetCell(0) );
				}
			}

			if ( pFace1->GetCellNum() > 1) {
				if ( pFace1->GetCell(0) == pCell ) {
					pCell1 = CAST(CCell*)( pFace1->GetCell(1) );
				} else {
					pCell1 = CAST(CCell*)( pFace1->GetCell(0) );
				}
			}

/*			if ( pCell0 != NULL) {
				if (pCell0->IsValid( pVert0, pFace->GetCollapsePos() ) == false )
					return lFIter;
			}

			if ( pCell1 != NULL) {
				if (pCell1->IsValid( pVert1, pFace->GetCollapsePos() ) == false )
					return lFIter;
			}
*/

			if ( pCell0 != NULL) {
				pCell0->RemoveCell( pCell );
				if ( pCell1 != NULL ) 
					pCell0->AddCell( pCell1 );
			}

			if ( pCell1 != NULL) {
				pCell1->RemoveCell( pCell );
				pCell1->RemoveFace( pFace1 );
				pCell1->AddFace( pFace0 );
				if ( pCell0 != NULL ) 
					pCell1->AddCell( pCell0 );
			}
			
			pFace0->RemoveCell( pCell );
			if ( pCell1 != NULL )
				pFace0->AddCell( pCell1 );
			pVert0->RemoveCell( pCell );
			pVert1->RemoveCell( pCell );
			pVert1->RemoveFace( pFace1 );
			if ( pVert1 == CAST(CVertex*)( pFace1->GetVert(0) ) ) {
				pVert1->RemoveVert( CAST(CVertex*)( pFace1->GetVert(1) ) );
				CAST(CVertex*)( pFace1->GetVert(1) )->RemoveVert( pVert1 );
				CAST(CVertex*)( pFace1->GetVert(1) )->RemoveFace( pFace1 );
				CAST(CVertex*)( pFace1->GetVert(1) )->RemoveCell( pCell );
			} else {
				pVert1->RemoveVert( CAST(CVertex*)( pFace1->GetVert(0) ) );
				CAST(CVertex*)( pFace1->GetVert(0) )->RemoveVert( pVert1 );
				CAST(CVertex*)( pFace1->GetVert(0) )->RemoveFace( pFace1 );
				CAST(CVertex*)( pFace1->GetVert(0) )->RemoveCell( pCell );
			}

			m_Faces.remove( pFace1 );
			delete pFace1;
			pFace1 = NULL;

			m_Cells.remove( pCell );
			delete pCell;
			pCell = NULL;
		}
		else { 	}  // more than 3 sides like quads
	}
	pVert1->RemoveVert( pVert0 );
	pVert1->RemoveFace( pFace );
	pVert0->RemoveVert( pVert1 );
	pVert0->RemoveFace( pFace );

	pVert0->SetPos( pFace->GetCollapsePos().GetPos() );

	for (i=0; i<pVert1->GetVertNum(); i++) {
		pVert0->AddVert( CAST(CVertex*)( pVert1->GetVert(i) ) );
		CAST(CVertex*)( pVert1->GetVert(i) )->RemoveVert( pVert1 );
		CAST(CVertex*)( pVert1->GetVert(i) )->AddVert( pVert0 );
	}
	for (i=0; i<pVert1->GetCellNum(); i++) {
		pCell = CAST(CCell*)( pVert1->GetCell(i) );
		pVert0->AddCell( CAST(CCell*)( pVert1->GetCell(i) ) );
		CAST(CCell*)( pVert1->GetCell(i) )->RemoveVert( pVert1 );
		CAST(CCell*)( pVert1->GetCell(i) )->AddVert( pVert0 );
	}
	for (i=0; i<pVert1->GetFaceNum(); i++) {
		pVert0->AddFace( CAST(CFace*)( pVert1->GetFace(i) ) );
		CAST(CFace*)( pVert1->GetFace(i) )->RemoveVert( pVert1 );
		CAST(CFace*)( pVert1->GetFace(i) )->AddVert( pVert0 );
	}

	for (i=0; i<pVert0->GetCellNum(); i++) {
		CAST(CCell*)( pVert0->GetCell(i) )->UpdateGeom() ;
	}

	for (i=0; i<pVert0->GetFaceNum(); i++) {
		CAST(CFace*)( pVert0->GetFace(i) )->UpdateGeom() ;
	}

	pFace->ClearNbrAll();
	for (m_FIter = m_Faces.begin(); m_FIter != m_Faces.end(); m_FIter++) {
		if ( (*m_FIter) == pFace ) {
			delete (*m_FIter);
			lFIter = m_Faces.erase( m_FIter );
			break;
		}
	}


	pVert1->ClearNbrAll();
	m_Verts.remove( pVert1 );
	delete pVert1;


	return lFIter;
}


//-------------------------------------
// Swap Face						  
//-------------------------------------
int CMeshAdap::SwapFace( CFace* pFace)
{
	int i;

	// ignore boundary faces
	if ( pFace->GetCellNum() < 2 ) 
		return 1;

	// only applicable to triangular cells
	if ( ( CAST(CCell*)( pFace->GetCell( 0 ) )->GetFaceNum() != 3 ) ||
		 ( CAST(CCell*)( pFace->GetCell( 1 ) )->GetFaceNum() != 3 ) ) 
		return 2;

	// create temporary cells
	int liNC =	2;
	vector<CCell*> lTmpCell( liNC );
	lTmpCell[0] = CAST(CCell*)( pFace->GetCell( 0 ) );
	lTmpCell[1] = CAST(CCell*)( pFace->GetCell( 1 ) );


	// cavity boundary vertices
	int liNV = 4;
	vector<CVertex*> lTmpBndV( liNV );	
	lTmpBndV[0] = CAST(CVertex*)( pFace->GetVert( 0 ) );
	lTmpBndV[1] = CAST(CVertex*)( pFace->GetVert( 1 ) );
	for ( i=0; i<3; i++) {
		if ( ( lTmpCell[0]->GetVert(i) != lTmpBndV[0] ) &&
			 ( lTmpCell[0]->GetVert(i) != lTmpBndV[1] ) )
			lTmpBndV[2] = CAST(CVertex*)( lTmpCell[0]->GetVert(i) );
		if ( ( lTmpCell[1]->GetVert(i) != lTmpBndV[0] ) &&
			 ( lTmpCell[1]->GetVert(i) != lTmpBndV[1] ) )
			lTmpBndV[3] = CAST(CVertex*)( lTmpCell[1]->GetVert(i) );
	}
	
	CFace* pF;
	// cavity boundary faces
	int liNF = 4;
	vector<CFace*> lTmpBndF( liNF );	
	for ( i=0; i<3; i++) {
		if  ( lTmpCell[0]->GetFace(i) != pFace ) {
			pF = CAST(CFace*)(lTmpCell[0]->GetFace(i));
			if ( CAST(CFace*)(lTmpCell[0]->GetFace(i))->IsNbrVert( lTmpBndV[0] ) )
				lTmpBndF[0] =  CAST(CFace*)( lTmpCell[0]->GetFace(i) );
			else
				lTmpBndF[1] =  CAST(CFace*)( lTmpCell[0]->GetFace(i) );
		}

		if  ( lTmpCell[1]->GetFace(i) != pFace ) {
			if ( CAST(CFace*)(lTmpCell[1]->GetFace(i))->IsNbrVert( lTmpBndV[0] ) )
				lTmpBndF[2] =  CAST(CFace*)( lTmpCell[1]->GetFace(i) );
			else
				lTmpBndF[3] =  CAST(CFace*)( lTmpCell[1]->GetFace(i) );
		}
	}

	// ignore concave quadrilaterals
	CVector2D lVec0 = lTmpBndV[1]->GetPoint().GetPos() - lTmpBndV[0]->GetPoint().GetPos();
	CVector2D lVec1 = lTmpBndV[2]->GetPoint().GetPos() - lTmpBndV[0]->GetPoint().GetPos();
	CVector2D lVec2 = lTmpBndV[3]->GetPoint().GetPos() - lTmpBndV[0]->GetPoint().GetPos();
	double ldTheta1 = lVec0.get_angle( lVec1 ) + lVec0.get_angle( lVec2 );
	if ( ldTheta1 > PI )
		return 3;
	// Delaunay
	ldTheta1 =  acos( ( ( pFace->GetMetric() * lVec1 ) * lVec2 ) /
					( sqrt( ( pFace->GetMetric() * lVec1 ) * lVec1 ) *
						sqrt( ( pFace->GetMetric() * lVec2 ) * lVec2 ) ) );


	lVec0 = lTmpBndV[0]->GetPoint().GetPos() - lTmpBndV[1]->GetPoint().GetPos();
	lVec1 = lTmpBndV[2]->GetPoint().GetPos() - lTmpBndV[1]->GetPoint().GetPos();
	lVec2 = lTmpBndV[3]->GetPoint().GetPos() - lTmpBndV[1]->GetPoint().GetPos();
	double ldTheta2 = lVec0.get_angle( lVec1 ) + lVec0.get_angle( lVec2 );
	if ( ldTheta2 > PI )
		return 3;
	// Delaunay
	ldTheta2 =  acos( ( ( pFace->GetMetric() * lVec1 ) * lVec2 ) /
				( sqrt( ( pFace->GetMetric() * lVec1 ) * lVec1 ) *
				  sqrt( ( pFace->GetMetric() * lVec2 ) * lVec2 ) ) );

	if ( (ldTheta1 + ldTheta2) > PI )
		return 4;
		
	pFace->ClearNbrVerts();
	pFace->AddVert( lTmpBndV[2] );
	pFace->AddVert( lTmpBndV[3] );

	lTmpBndV[0]->RemoveCell( lTmpCell[1] );
	lTmpBndV[0]->RemoveFace( pFace );
	lTmpBndV[0]->RemoveVert( lTmpBndV[1] );

	lTmpBndV[1]->RemoveCell( lTmpCell[0] );
	lTmpBndV[1]->RemoveFace( pFace );
	lTmpBndV[1]->RemoveVert( lTmpBndV[0] );

	lTmpBndV[2]->AddCell( lTmpCell[1] );
	lTmpBndV[2]->AddFace( pFace );
	lTmpBndV[2]->AddVert( lTmpBndV[3] );
	
	lTmpBndV[3]->AddCell( lTmpCell[0] );
	lTmpBndV[3]->AddFace( pFace );
	lTmpBndV[3]->AddVert( lTmpBndV[2] );

	lTmpCell[0]->RemoveVert( lTmpBndV[1] );
	lTmpCell[0]->AddVert( lTmpBndV[3] );
	
	lTmpCell[1]->RemoveVert( lTmpBndV[0] );
	lTmpCell[1]->AddVert( lTmpBndV[2] );
	
	lTmpCell[0]->RemoveFace( lTmpBndF[1] );
	lTmpCell[0]->AddFace( lTmpBndF[2] );
	
	lTmpCell[1]->RemoveFace( lTmpBndF[2] );
	lTmpCell[1]->AddFace( lTmpBndF[1] );
	
	if ( lTmpBndF[2]->GetCellNum() > 1 ) {
		if ( lTmpBndF[2]->GetCell( 0 ) == lTmpCell[1] ){
			CAST(CCell*)( lTmpBndF[2]->GetCell(1) )->RemoveCell( lTmpCell[1] );
			CAST(CCell*)( lTmpBndF[2]->GetCell(1) )->AddCell( lTmpCell[0] );
			lTmpCell[1]->RemoveCell( lTmpBndF[2]->GetCell(1) );
			lTmpCell[0]->AddCell( lTmpBndF[2]->GetCell(1) );
		} else {
			CAST(CCell*)( lTmpBndF[2]->GetCell(0) )->RemoveCell( lTmpCell[1] );
			CAST(CCell*)( lTmpBndF[2]->GetCell(0) )->AddCell( lTmpCell[0] );
			lTmpCell[1]->RemoveCell( lTmpBndF[2]->GetCell(0) );
			lTmpCell[0]->AddCell( lTmpBndF[2]->GetCell(0) );
		}
	}

	if ( lTmpBndF[1]->GetCellNum() > 1 ) {
		if ( lTmpBndF[1]->GetCell( 0 ) == lTmpCell[0] ){
			CAST(CCell*)( lTmpBndF[1]->GetCell(1) )->RemoveCell( lTmpCell[0] );
			CAST(CCell*)( lTmpBndF[1]->GetCell(1) )->AddCell( lTmpCell[1] );
			lTmpCell[0]->RemoveCell( lTmpBndF[1]->GetCell(1) );
			lTmpCell[1]->AddCell( lTmpBndF[1]->GetCell(1) );
		} else {
			CAST(CCell*)( lTmpBndF[1]->GetCell(0) )->RemoveCell( lTmpCell[0] );
			CAST(CCell*)( lTmpBndF[1]->GetCell(0) )->AddCell( lTmpCell[1] );
			lTmpCell[0]->RemoveCell( lTmpBndF[1]->GetCell(0) );
			lTmpCell[1]->AddCell( lTmpBndF[1]->GetCell(0) );
		}
	}

	lTmpBndF[1]->RemoveCell( lTmpCell[0] );
	lTmpBndF[1]->AddCell( lTmpCell[1] ); 
	
	lTmpBndF[2]->RemoveCell( lTmpCell[1] );
	lTmpBndF[2]->AddCell( lTmpCell[0] ); 

	for (i=0; i<liNC; i++)
		lTmpCell[i]->UpdateGeom();

	for (i=0; i<liNF; i++)
		lTmpBndF[i]->UpdateGeom();

	pFace->UpdateGeom();

	m_iSwapFace++;
	return 0;
}


void CMeshAdap::Refine()
{
	int n ;
	double ldSize;
	CVector2D lVec;

	for ( n = 0; n != m_Faces.size();  ) 
	{
		n = m_Faces.size();
		for ( m_FIter = m_Faces.begin(); m_FIter != m_Faces.end(); m_FIter++) {
			lVec = (*m_FIter)->GetUTV() * (*m_FIter)->GetArea();
			ldSize = ( (*m_FIter)->GetMetric() * lVec ) * lVec;
			if ( ldSize > 2.0 )
				SplitFace( *m_FIter );
		}
	}
}

void CMeshAdap::Coarsen()
{
	int n; 
	double ldSize;
	CVector2D lVec;

	for ( n = 0 ; n != m_Faces.size(); ) 
	{
		n = m_Faces.size();
		for ( m_FIter = m_Faces.begin(); m_FIter != m_Faces.end();) 	{
			lVec = (*m_FIter)->GetUTV() * (*m_FIter)->GetArea();
			ldSize = ( (*m_FIter)->GetMetric() * lVec ) * lVec;
			if ( ldSize < 0.5 ) 
				m_FIter = CollapseFace( m_FIter );
			else
				m_FIter++;
		}
	}
}


void CMeshAdap::Swap()
{
	m_iSwapFace = 1;
	for ( ; m_iSwapFace != 0; ) {
		m_iSwapFace = 0;
		for ( m_FIter = m_Faces.begin(); m_FIter != m_Faces.end(); m_FIter++) 
			SwapFace( *m_FIter );
	}
}


void CMeshAdap::Smooth()
{
	int i;
	int liNV;	// number of neighbor vertices
	CVector2D lVec;


	for (m_VIter = m_Verts.begin(); m_VIter != m_Verts.end(); m_VIter++) {
		// interior points
		if ( (*m_VIter)->GetType() == 0) {
		liNV = (*m_VIter)->GetVertNum();
		lVec.SetVector( 0.0, 0.0 );
		for ( i=0; i<liNV; i++) {
			lVec = lVec + static_cast<CVertex*>( (*m_VIter)->GetVert( i ) )->GetPoint().GetPos();
		}	
		lVec /= static_cast<double>( liNV );
		lVec = lVec *0.5 + (*m_VIter)->GetPoint().GetPos() * 0.5;
		(*m_VIter)->SetPos( lVec );
		} 
	}
}


void CMeshAdap::SmoothAnis()
{
	int i;
	int liNF;	// number of neighbor vertices
	CVector2D lVec_i;
	CVector2D lVec;
	double ldK, ldK_i;
	CFace* lpFace;


	for (m_VIter = m_Verts.begin(); m_VIter != m_Verts.end(); m_VIter++) {
		// interior points
		if ( (*m_VIter)->GetType() == 0) {
		liNF = (*m_VIter)->GetFaceNum();
		lVec.SetVector( 0.0, 0.0 );
		ldK = 0.0;
		for ( i=0; i<liNF; i++) {
			lpFace = CAST(CFace*)( (*m_VIter)->GetFace( i ) );
			if ( lpFace->GetVert(0) == (*m_VIter) ) 
				lVec_i = CAST(CVertex*)(lpFace->GetVert(1))->GetPoint().GetPos();
			else
				lVec_i = CAST(CVertex*)(lpFace->GetVert(0))->GetPoint().GetPos();
			ldK_i = lpFace->GetTargetArea() / lpFace->GetArea();

			lVec += lVec_i * ldK_i;
			ldK += ldK_i; 
		}	
		lVec /= ldK;
		lVec = lVec *1.0 + (*m_VIter)->GetPoint().GetPos() * 0.0;
		(*m_VIter)->SetPos( lVec );
		} 
	}
}

void CMeshAdap::Topolog()
{
	int n;
	CVertex* lpVert0;
	CVertex* lpVert1;

	n = 0;
	for ( ; n != m_Faces.size(); ) {
		n = m_Faces.size();
		for ( m_FIter = m_Faces.begin(); m_FIter != m_Faces.end();) 	{
			lpVert0 = CAST(CVertex*)( (*m_FIter)->GetVert(0) );
			lpVert1 = CAST(CVertex*)( (*m_FIter)->GetVert(1) );

			if ( ( (lpVert0->GetFaceNum() + lpVert1->GetFaceNum()) == 10) &&
				 (lpVert0->GetType() ==0 ) && (lpVert1->GetType() ==0 ) )
				m_FIter = CollapseFace( m_FIter );
			else
				m_FIter++;
		}
	}
}
