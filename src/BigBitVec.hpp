/*
 * Copyright (C) 2006-2007 Chris Hamilton <chamilton@cs.dal.ca>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef _BIGBITVEC_HPP_
#define _BIGBITVEC_HPP_

#include "Rcpp.h"
#include "Common.hpp"
#include "FixBitVec.hpp"
#include <string.h>


#define BBV_MIN(a,b)				((a)<(b)?(a):(b))
#define BBV_MAX(a,b)				((a)>(b)?(a):(b))
#define FBVS_NEEDED(b)			((BBV_MAX(b,1)+FBV_BITS-1)/FBV_BITS)
#define BBV_MODSPLIT(r,b,k) { b=(k); r=b/FBV_BITS; b-=r*FBV_BITS; }


class CBigBitVec
{
	public:

		static
		EBitVecType
		Type()
		{
			return eBig;
		}

		// Constructor, with optional number of bits.
		H_INLINE
		CBigBitVec(
			int iBits = FBV_BITS
			)
		{
			// Determine number of racks required.
			m_iRacks = FBVS_NEEDED(iBits);
		
			// Allocate the memory.
			m_pcRacks = new CFixBitVec[m_iRacks];
					
			return;
		}
				
				
		// Copy construct.  Creates duplicate.
		H_INLINE
		CBigBitVec(
			const CBigBitVec &cBBV
			)
		{
			int i;
		
			m_iRacks = cBBV.m_iRacks;
			m_pcRacks = new CFixBitVec[m_iRacks];
		
			// Copy the rack values.
			/*for ( i = 0; i < m_iRacks; i++ )
				m_pcRacks[i] = cBBV.m_pcRacks[i];*/
			memcpy( static_cast<void*>(m_pcRacks),
				static_cast<const void*>(cBBV.m_pcRacks),
				sizeof(CFixBitVec)*m_iRacks );
		
			return;
		}
		
		
		// Copy constructor.
		H_INLINE
		CBigBitVec(
			const CFixBitVec &cFBV
			)
		{
			m_iRacks = 1;
			m_pcRacks = new CFixBitVec[m_iRacks];
		
			m_pcRacks[0] = cFBV;
		
			return;
		}
		
		
		// Destructor
		H_INLINE
		~CBigBitVec()
		{
			delete [] m_pcRacks;
		
			return;
		}
		
		
		// Returns the current size in bits.
		H_INLINE
		int
		GetSize() const
		{
			return m_iRacks*FBV_BITS;
		}
		
		
		// Resize function.  Returns the number of bits
		// we can accomodate after resizing.
		H_INLINE
		CBigBitVec &
		SetSize(
			int iBits
			)
		{
			int i;
			
			// How many racks do we need?
			int iRacks = FBVS_NEEDED(iBits);
		
			// Same size?
			if ( iRacks == m_iRacks )
				return (*this);
			
			// Allocate new racks.
			CFixBitVec *pcRacks = new CFixBitVec[iRacks];
		
			// Copy over the old values.
			/*for ( i = 0; i < BBV_MIN(iRacks,m_iRacks); i++ )
				pcRacks[i] = m_pcRacks[i];*/
			i = BBV_MIN(iRacks,m_iRacks);
			memcpy( static_cast<void*>(pcRacks),
				static_cast<const void*>(m_pcRacks),
				sizeof(CFixBitVec)*i );
		
			// Zero the new values
			/*for ( ; i < iRacks; i++ )
				pcRacks[i].Zero();*/
			if ( iRacks > i )
				memset( static_cast<void*>(pcRacks + i), 0,
					sizeof(CFixBitVec)*(iRacks-i) );
		
			// Release the old stuff.
			delete [] m_pcRacks;
		
			// Replace old with new.
			m_iRacks = iRacks;
			m_pcRacks = pcRacks;
			
			return (*this);
		}
		
		
		// Zeros the bit-vector.
		H_INLINE
		CBigBitVec &
		Zero()
		{
			int i;
			
			/*for ( i = 0; i < m_iRacks; i++ )
				m_pcRacks[i].Zero();*/
			memset( static_cast<void*>(m_pcRacks), 0,
				sizeof(CFixBitVec)*m_iRacks );
		
			return (*this);
		}
		
		
		// Truncates the bit-vector to a given precision in
		// bits (zeroes MSBs without shrinking the vector)
		H_INLINE
		CBigBitVec &
		Truncate(
			int iBits
			)
		{
			assert( iBits >= 0 && iBits <= GetSize() );
			int r, b, i;
		
			BBV_MODSPLIT(r,b,iBits);
			
			if ( r >= m_iRacks )
				return (*this);
		
			m_pcRacks[r].Truncate(b);
			
			for ( i = r+1; i < m_iRacks; i++ )
				m_pcRacks[i].Zero();
		
			return (*this);
		}
		
		
		// Assignment operator.  No resizing.
		H_INLINE
		CBigBitVec &
		operator=(
			const CBigBitVec &cBBV
			)
		{
			int i;
			
			if ( m_iRacks < cBBV.m_iRacks )
			{
				/*for ( i = 0; i < m_iRacks; i++ )
					m_pcRacks[i] = cBBV.m_pcRacks[i];*/
				memcpy( static_cast<void*>(m_pcRacks),
					static_cast<const void*>(cBBV.m_pcRacks),
					sizeof(CFixBitVec)*m_iRacks );
			}
			else
			{
				/*for ( i = 0; i < cBBV.m_iRacks; i++ )
					m_pcRacks[i] = cBBV.m_pcRacks[i];
				for ( ; i < m_iRacks; i++ )
					m_pcRacks[i].Zero();*/
				memcpy( static_cast<void*>(m_pcRacks),
					static_cast<const void*>(cBBV.m_pcRacks),
					sizeof(CFixBitVec)*cBBV.m_iRacks );
				memset( static_cast<void*>(m_pcRacks+cBBV.m_iRacks),
					0, sizeof(CFixBitVec)*(m_iRacks-cBBV.m_iRacks) );
			}
		
			return (*this);
		}
		H_INLINE
		CBigBitVec &
		operator=(
			const CFixBitVec &cFBV
			)
		{
			int i;
			m_pcRacks[0] = cFBV;
			/*for ( i = 1; i < m_iRacks; i++ )
				m_pcRacks[i].Zero();*/
			memset( static_cast<void*>(m_pcRacks+1),
				0, sizeof(CFixBitVec)*(m_iRacks-1) );
			return (*this);
		}
		H_INLINE
		CBigBitVec &
		operator=(
			 FBV_UINT j
			)
		{
			int i;
			m_pcRacks[0] = j;
			/*for ( i = 1; i < m_iRacks; i++ )
				m_pcRacks[i].Zero();*/
			memset( static_cast<void*>(m_pcRacks+1),
				0, sizeof(CFixBitVec)*(m_iRacks-1) );
			return (*this);
		}
		
		// Returns the value of the nth bit.
		H_INLINE
		bool
		GetBit(
			int iIndex
			) const
		{
			assert( iIndex >= 0 && iIndex < GetSize() );
			int r, b;
			BBV_MODSPLIT(r,b,iIndex);
			return m_pcRacks[r].GetBit(b);
		}
		
		
		// Sets the value of the nth bit.
		H_INLINE
		CBigBitVec &
		SetBit(
			int iIndex,
			bool bBit
			)
		{
			assert( iIndex >= 0 && iIndex < GetSize() );
			int r, b;
			BBV_MODSPLIT(r,b,iIndex);
			m_pcRacks[r].SetBit(b,bBit);
			return (*this);
		}
		
		
		// Toggles the value of the nth bit.
		H_INLINE
		CBigBitVec &
		ToggleBit(
			int iIndex
			)
		{
			assert( iIndex >= 0 && iIndex < GetSize() );
			int r, b;
			BBV_MODSPLIT(r,b,iIndex);
			m_pcRacks[r].ToggleBit(b);
			return (*this);
		}
		
		
		// In place AND.
		H_INLINE
		CBigBitVec &
		operator&=(
			const CBigBitVec &cBBV
			)
		{
			int i;
			
			for ( i = 0; i < BBV_MIN(m_iRacks,cBBV.m_iRacks); i++ )
				m_pcRacks[i] &= cBBV.m_pcRacks[i];
			
			return (*this);
		}
		H_INLINE
		CBigBitVec &
		operator&=(
			const CFixBitVec &r
			)
		{
			m_pcRacks[0] &= r;
			return (*this);
		}
		H_INLINE
		CBigBitVec &
		operator&=(
			 FBV_UINT i
			)
		{
			m_pcRacks[0] &= i;
			return (*this);
		}
		
		
		// AND operator.
		H_INLINE
		CBigBitVec
		operator&(
			const CBigBitVec &cBBV
			) const
		{
			CBigBitVec t( *this );
			t &= cBBV;
		
			return t;
		}
		H_INLINE
		CBigBitVec
		operator&(
			const CFixBitVec &r
			)
		{
			CBigBitVec t( *this );
			t &= r;
			
			return t;
		}
		H_INLINE
		CBigBitVec
		operator&(
			 FBV_UINT i
			)
		{
			CBigBitVec t( *this );
			t &= i;
			
			return t;
		}
		
		
		// In place OR.
		H_INLINE
		CBigBitVec &
		operator|=(
			const CBigBitVec &cBBV
			)
		{
			int i;
			
			for ( i = 0; i < BBV_MIN(m_iRacks,cBBV.m_iRacks); i++ )
				m_pcRacks[i] |= cBBV.m_pcRacks[i];
			
			return (*this);
		}
		H_INLINE
		CBigBitVec &
		operator|=(
			const CFixBitVec &r
			)
		{
			m_pcRacks[0] |= r;
			return (*this);
		}
		H_INLINE
		CBigBitVec &
		operator|=(
			 FBV_UINT i
			)
		{
			m_pcRacks[0] |= i;
			return (*this);
		}
		
		
		// OR operator.
		H_INLINE
		CBigBitVec
		operator|(
			const CBigBitVec &cBBV
			) const
		{
			CBigBitVec t( *this );
			t |= cBBV;
		
			return t;
		}
		H_INLINE
		CBigBitVec
		operator|(
			const CFixBitVec &r
			)
		{
			CBigBitVec t( *this );
			t |= r;
			
			return t;
		}
		H_INLINE
		CBigBitVec
		operator|(
			 FBV_UINT i
			)
		{
			CBigBitVec t( *this );
			t |= i;
			
			return t;
		}
		
		
		// In place XOR.
		H_INLINE
		CBigBitVec &
		operator^=(
			const CBigBitVec &cBBV
			)
		{
			int i;
			
			for ( i = 0; i < BBV_MIN(m_iRacks,cBBV.m_iRacks); i++ )
				m_pcRacks[i] ^= cBBV.m_pcRacks[i];
			
			return (*this);
		}
		H_INLINE
		CBigBitVec &
		operator^=(
			const CFixBitVec &r
			)
		{
			m_pcRacks[0] ^= r;
			return (*this);
		}
		H_INLINE
		CBigBitVec &
		operator^=(
			 FBV_UINT i
			)
		{
			m_pcRacks[0] ^= i;
			return (*this);
		}
		
		
		// XOR operator.
		H_INLINE
		CBigBitVec
		operator^(
			const CBigBitVec &cBBV
			) const
		{
			CBigBitVec t( *this );
			t ^= cBBV;
		
			return t;
		}
		H_INLINE
		CBigBitVec	
		operator^(
			const CFixBitVec &r
			)
		{
			CBigBitVec t( *this );
			t ^= r;
			
			return t;
		}
		H_INLINE
		CBigBitVec
		operator^(
			 FBV_UINT i
			)
		{
			CBigBitVec t( *this );
			t ^= i;
			
			return t;
		}
		
		
		// Shift left operation, in place.
		H_INLINE
		CBigBitVec &
		operator<<=(
			int iBits
			)
		{
			assert( iBits >= 0 );
			int r, b, i;
		
			// No shift?
			if ( iBits == 0 )
				return (*this);
		
			BBV_MODSPLIT(r,b,iBits);
		
			// All racks?
			if ( r >= m_iRacks )
			{
				Zero();
				return (*this);
			}
		
			// Do rack shifts.
			if ( r > 0 )
			{
				for ( i = m_iRacks-1; i >= r; i-- )
					m_pcRacks[i] = m_pcRacks[i-r];
				for ( ; i >= 0; i-- )
					m_pcRacks[i].Zero();
			}
		
			// Do bit shifts.
			if ( b > 0 )
			{
				int bi = FBV_BITS - b;
		    for ( i = m_iRacks-1; i >= r+1; i-- )
		    {
			    m_pcRacks[i] <<= b;
					m_pcRacks[i] |= m_pcRacks[i-1] >> bi;
			  }
			  m_pcRacks[i] <<= b;
			}
		
			return (*this);
		}
		
		
		// Shift left operation.
		H_INLINE
		CBigBitVec
		operator<<(
			int iBits
			) const
		{
			CBigBitVec t( *this );
			t <<= iBits;
			return t;
		}
		
		
		// Shift right operation, in place.
		H_INLINE
		CBigBitVec &
		operator>>=(
			int iBits
			)
		{
			assert( iBits >= 0 );
			int r, b, i;
		
			// No shift?
			if ( iBits == 0 )
				return (*this);
		
			BBV_MODSPLIT(r,b,iBits);
		
			// All racks?
			if ( r >= m_iRacks )
			{
				Zero();
				return (*this);
			}
		
			// Do rack shifts.
			if ( r > 0 )
			{
				for ( i = 0; i < m_iRacks-r; i++ )
					m_pcRacks[i] = m_pcRacks[i+r];
				for ( ; i < m_iRacks; i++ )
					m_pcRacks[i].Zero();
			}
		
			// Do bit shifts.
			if ( b > 0 )
			{
				int bi = FBV_BITS - b;
		    for ( i = 0; i < m_iRacks-r-1; i++ )
		    {
			    m_pcRacks[i] >>= b;
					m_pcRacks[i] |= m_pcRacks[i+1] << bi;
			  }
			  m_pcRacks[i] >>= b;
			}
		
			return (*this);
		}
		
		
		// Shift right operation.
		H_INLINE
		CBigBitVec
		operator>>(
			int iBits
			) const
		{
			CBigBitVec t( *this );
			t >>= iBits;
			return t;
		}
		
		
		// Right rotation, in place.
		H_INLINE
		CBigBitVec &
		Rotr(
			int iBits,
			int iWidth = 0
			)
		{
			assert( iBits >= 0 );
			int r, bl, b, i;
		
			// Fill in the width, if necessary.
			if ( iWidth <= 0 )
				iWidth = GetSize();
		
			// Modulo the number of bits.
			//#D FBVMOD(iBits,iWidth);
			assert( iBits < iWidth );
			if ( iBits == 0 ) return (*this);
		
			// Ensure we are truncated appropriately.
			Truncate(iWidth);
		
			CBigBitVec t1( *this );
			(*this) >>= iBits;
			t1 <<= (iWidth-iBits);
			(*this) |= t1;
			
			Truncate(iWidth);
		
			return (*this);
		}
		
		
		// Right rotation.
		H_INLINE
		CBigBitVec
		RotrCopy(
			int iBits,
			int iWidth = 0
			) const
		{
			CBigBitVec t( *this );
			t.Rotr(iBits,iWidth);
			return t;
		}
		
		
		// Left rotation, in place.
		H_INLINE
		CBigBitVec &
		Rotl(
			int iBits,
			int iWidth = 0
			)
		{
			assert( iBits >= 0 );
			int r, bl, b, i;
		
			// Fill in the width, if necessary.
			if ( iWidth <= 0 )
				iWidth = GetSize();
		
			// Modulo the number of bits.
			//FBVMOD(iBits,iWidth);
			assert( iBits < iWidth );
			if ( iBits == 0 ) return (*this);
		
			// Ensure we are truncated appropriately.
			Truncate(iWidth);
		
			CBigBitVec t1( *this );
			(*this) <<= iBits;
			t1 >>= (iWidth-iBits);
			(*this) |= t1;

			Truncate(iWidth);
		
			return (*this);
		}
		
		
		// Left rotation.
		H_INLINE
		CBigBitVec
		RotlCopy(
			int iBits,
			int iWidth = 0
			) const
		{
			CBigBitVec t( *this );
			t.Rotl(iBits,iWidth);
			return t;
		}
		
		
		// Returns true if the rack is zero valued.
		H_INLINE
		bool
		IsZero() const
		{
			int i;
			for ( i = 0; i < m_iRacks; i++ )
				if ( !m_pcRacks[i].IsZero() ) return false;
			return true;
		}
		
		
		// Returns the number of trailing set bits.
		H_INLINE
		int
		Tsb() const
		{
			int c, i, j;
			c = 0;
			for ( i = 0; i < m_iRacks; i++ )
			{
				j = m_pcRacks[i].Tsb();
				c += j;
				if ( j < FBV_BITS )
					break;
			}
			return c;
		}
		
		
		// Returns the index of the first set bit.
		// (numbered 1 to n, with 0 meaning no bits were set)
		H_INLINE
		int
		Fsb() const
		{
			int c, i, j;
			c = 0;
			for ( i = 0; i < m_iRacks; i++ )
			{
				j = m_pcRacks[i].Fsb();
				if ( j < FBV_BITS )
					return c + j;
				else
					c += FBV_BITS;
			}
			return 0;
		}
		
		
		// Prefix decrement.  Returns true if there
		// was a carry, false otherwise.
		H_INLINE
		bool
		operator--()
		{
			int i = 0;
			bool b;
			while ( i < m_iRacks && (b = --m_pcRacks[i]) ) i++;
			
			return b;
		}
		
		
		// Gray Code
		H_INLINE
		CBigBitVec &
		GrayCode()
		{
			int i;
			FBV_UINT s = 0;
		
			for ( i = m_iRacks-1; i >= 0; i-- )
			{
				FBV_UINT t = m_pcRacks[i].Rack() & 1;
				m_pcRacks[i].GrayCode();
				m_pcRacks[i] ^= (s<<(FBV_BITS-1));
				s = t;
			}
		
			return (*this);
		}
		
		
		// Gray Code Inverse
		H_INLINE
		CBigBitVec &
		GrayCodeInv()
		{
			int i;
			FBV_UINT s = 0;
		
			for ( i = m_iRacks-1; i >= 0; i-- )
			{
				m_pcRacks[i].GrayCodeInv();
				if ( s ) m_pcRacks[i].Complement();
				s = m_pcRacks[i].Rack() & 1;
			}
		}
		
		
		// Complement
		H_INLINE
		CBigBitVec &
		Complement()
		{
			int i;
			for ( i = 0; i < m_iRacks; i++ )
				m_pcRacks[i].Complement();
			return (*this);
		}
		
		
		// Returns the first rack.
		H_INLINE
		FBV_UINT &
		Rack()
		{
			return m_pcRacks[0].Rack();
		}
		H_INLINE
		FBV_UINT
		Rack() const
		{
			return m_pcRacks[0].Rack();
		}
		
		
		// Returns the racks.
		H_INLINE
		FBV_UINT *
		Racks()
		{
			return reinterpret_cast<FBV_UINT*>(m_pcRacks);
		}
		H_INLINE
		const FBV_UINT *
		Racks() const
		{
			return reinterpret_cast<FBV_UINT*>(m_pcRacks);
		}
		
		
		// Returns the number of racks.
		H_INLINE
		int
		RackCount() const
		{
			return m_iRacks;
		}

	
	private:	
		

		// Right rotates entire racks (in place).
		H_INLINE
		void
		RackRotr(
			int k
			)
		{
			assert( 0 <= k && k < m_iRacks );
			
			int c, v;
			CFixBitVec tmp;
		
		  if (k == 0) return;
		
		  c = 0;
		  for (v = 0; c < m_iRacks; v++)
			{
		    int t = v, tp = v + k;
		    tmp = m_pcRacks[v];
		    c++;
		    while (tp != v)
				{
		      m_pcRacks[t] = m_pcRacks[tp];
		      t = tp;
		      tp += k;
		      if (tp >= m_iRacks) tp -= m_iRacks;
		      c++;
		    }
		    m_pcRacks[t] = tmp;
		  }
		
		  return;
		}


		CFixBitVec *m_pcRacks;
		int m_iRacks;
};


#endif
