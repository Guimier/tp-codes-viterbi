#include <iostream>       
#include <fstream>
#include <string>        
#include <vector>
#include <bitset>        
#include <cstdlib>
#include <ctime>
#include <limits>

const double errorRate = 0.1;

const int N=2;
const int K=1;
const int R=4;

#undef DEBUG

using namespace std; 

////////////////////////////////////////////////////////////
//      template<int bits> bitset<bits> randBitset()      //
//                                                        //
//               Generate random bitset                   //
////////////////////////////////////////////////////////////

template<int bits>
bitset<bits> randBitset()
{
	bitset<bits> r( rand() );
	for ( int i = 0; i < bits/16 - 1; i++ ) {
		r <<= 16;
		r |= bitset<bits>( rand() );
	}
	return r;
}

////////////////////////////////////////////////////////////
// vector< bitset<N> > GSM_code(vector< bitset<K> > mess) //
//                                                        //
//     Convolutional coding of a message (GSM norm)       //
////////////////////////////////////////////////////////////

vector< bitset<N> > GSM_code( vector< bitset<K> > mess )
{
	int g0, g1;
	vector< bitset<N> > mess_out;

	bitset<N> cod_out; 
	bitset<R+1> G0(25);
	bitset<R+1> G1(27); 
	bitset<R+1> reg; 
	reg.reset();
 
	#ifdef DEBUG
		cout << "-------------------- Debug Informations (Coding) --------------------" << endl << endl;
		cout << "Initial register ( u(i-4)  u(i-3)  u(i-2)  u(i-1)  u(i)  ): " << reg << endl;
		cout << "Polynom G0       ( g0(i-4) g0(i-3) g0(i-2) g0(i-1) g0(i) ): " << G0 << endl;
		cout << "Polynom G1       ( g1(i-4) g1(i-3) g1(i-2) g1(i-1) g1(i) ): " << G1 << endl << endl;
	#endif

	for ( vector<bitset<K> >::iterator it = mess.begin(); it != mess.end(); ++it ) {
		reg = reg << 1;
		reg.set( 0, (*it).count() );

		g0 = ( reg & G0 ).count() % 2;
		g1 = ( reg & G1 ).count() % 2;

		cod_out.reset();
		cod_out.set( 0, g0 );
		cod_out.set( 1, g1 );

		mess_out.push_back( cod_out );

		#ifdef DEBUG
			cout << "Block number: " << ++i << " - In frame: "<< *it << endl; 
			cout << "\t Current status of registers: "<< reg << endl;
			cout << "\t Out : " << cod_out << endl;
		#endif
	}
	
	#ifdef DEBUG
		cout << "------------------------------------------------------------------" << endl << endl;
	#endif

	return mess_out;
}

bitset<N> automatCycle( bool in, bitset<R> registers )
{
	bitset<N> cod_out( 0 );
	bitset<R+1> G0(25);
	bitset<R+1> G1(27);
	bitset<R+1> reg( registers.to_ulong() << 1 );
	
	reg.set( 0, in );

	int g0 = ( reg & G0 ).count() % 2;
	int g1 = ( reg & G1 ).count() % 2;

	cod_out.set( 0, g0 );
	cod_out.set( 1, g1 );

	return cod_out;
}

/////////////////////////////////////////////////////////////////////////
// vector< bitset<N> >  GSM_transmission(vector< bitset<N> > mess_cod) //
//                                                                     //
//         Simulation of a transmission channel => adding errors       //
/////////////////////////////////////////////////////////////////////////

vector< bitset<N> > GSM_transmission( vector< bitset<N> > mess_cod )
{
	vector< bitset<N> > mess_tra = mess_cod;

	int errCount = 0;
	
	for ( vector< bitset<N> >::iterator it = mess_tra.begin(); it != mess_tra.end(); ++it ) {
		for ( int i = 0 ; i < N; ++i ) {
			if ( rand() <= RAND_MAX * errorRate ) {
				++errCount;
				(*it)[i] = ! (*it)[i];
			}
		}
	}
	
	cout << "Error count      : " << errCount << endl;

	return mess_tra;
}

/** A path of bits. */
class Path: public vector<bool> {
	
	public:
	
		/** Create an empty path. */
		Path(): vector<bool>( 0 )
		{}
		
		/** Create a path going one step further.
		 * @param prev
		 *     Reference path.
		 * @param last
		 *     Additional bit.
		 */
		Path( const Path& prev, const bool last ):
			vector<bool>( prev )
		{
			push_back( last );
		}
		
		/** Copy a path.
		 * @param p
		 *     Path to be copied.
		 */
		Path( const Path& p ):
			vector<bool>( p )
		{}
	
};

/** Path with an associated weight. */
class WeightedPath {
	
	public:
		/** Integral type for weight */
		typedef int weight_t;
		/** Weight considered infinite. */
		static const weight_t INFINITY;
	
	private:
		/** Path. */
		Path path;
		/** Weight. */
		weight_t weight;
	
	public:
		/** Create weighted path.
		 * @param path
		 *     Path.
		 * @param weight
		 *     Weight.
		 */
		WeightedPath( const Path& path, const weight_t weight ):
			path( path ),
			weight( weight )
		{}
		/** Create weighted path going one step further.
		 * @param prev
		 *     Reference path.
		 * @param last
		 *     Additional bit.
		 * @param add
		 *     Additional weight.
		 */
		WeightedPath( const WeightedPath& prev, bool last, const weight_t add ):
			path( prev.path, last ),
			weight( INFINITY - prev.weight <= add ? INFINITY : prev.weight + add )
		{}
		/** Copy weighted path.
		 * @param wp 
		 *     Weighted path to be copied.
		 */
		WeightedPath( const WeightedPath& wp ):
			path( wp.path ),
			weight( wp.weight )
		{}
		/** Create an empty, no-weight path. */
		WeightedPath():
			path( Path() ),
			weight( 0 )
		{}
		
		/** Get the weight. */
		weight_t getWeight()
		{
			return weight;
		}
		
		/** Get the path. */
		Path getPath()
		{
			return path;
		}
	
};

const WeightedPath::weight_t WeightedPath::INFINITY = std::numeric_limits<weight_t>::max();

/** Compute the Hamming distance between two codes.
 * @param W
 *    Size of code we compare
 * @param a, b
 *	Codes to be compared
 */
template <size_t W>
size_t hammingDistance( const bitset<W>& a, const bitset<W>& b )
{
	return ( a ^ b ).count();
}

//////////////////////////////////////////////////////////////////
// vector< bitset<K> > GSM_decode(vector< bitset<N> > mess_tra) //
//                                                              //
//     Convolutional decoding of a message (GSM norm)           //
//////////////////////////////////////////////////////////////////
vector< bitset<K> > GSM_decode( vector< bitset<N> > mess_tra )
{
	vector< bitset<K> > mess_dec;
	
	// Template used for impossible paths initialization.
	const WeightedPath impossiblePath( Path(), WeightedPath::INFINITY );
	// Paths associated with the states *before* the new bits.
	vector<WeightedPath> prevPaths;
	// State ( 0, 0, 0, 0 ) is possible when the algorithm starts, 0 errors.
	prevPaths.push_back( WeightedPath( Path(), 0 ) );
	// All other states are impossible.
	for ( int i = 1; i < 1 << R; ++i ) {
		prevPaths.push_back( impossiblePath );
	}
	
	// Paths associated with the states *after* the new bits.
	vector<WeightedPath> nextPaths( 1 << R );
	
	// Precomputed theoric emissions of the automat.
	// First index: automate *output* state.
	// Second index: bit the automat have forgotten in the operation.
	vector< vector< bitset<N> > > theoricEmmissions;
	
	// Initialize matrix intermediate vectors.
	vector< bitset<N> > templ( 2 );
	for ( int i = 0; i < 1 << R; ++i ) {
		theoricEmmissions.push_back( templ );
	}
	
	// For all automat *output* states.
	for ( int i = 0; i < ( 1 << R ); ++i ) {
		// For all the bits *forgotten* by the automat.
		for ( int j = 0; j < ( 1 << K ); ++j ) {
			 theoricEmmissions[i][j] = automatCycle(
				 // Added bit
				i & 1,
				// *Input* state.
				bitset<R>( ( i >> 1 ) | ( j << ( R - 1 ) ) )
			 );
		}
	}
	
	// For all couples of received bits.
	for ( vector< bitset<N> >::const_iterator it = mess_tra.begin(); it != mess_tra.end(); ++it ) {
		// For all automat *output* states.
		for ( int i = 0; i < ( 1 << R ); ++i ) {
			WeightedPath forgotZero(
				// Best path in state before the automat forgets a zero.
				prevPaths[i>>1],
				// Bit added.
				i & 1,
				hammingDistance( theoricEmmissions[i][0], *it )
			);
			
			WeightedPath forgotOne(
				// Best path in state before the automat forgets a one.
				prevPaths[(i>>1)|(1<<(R-1))],
				// Bit added.
				i & 1,
				hammingDistance( theoricEmmissions[i][1], *it )
			);
			
			// Remember best path for this *output* state.
			if ( forgotZero.getWeight() < forgotOne.getWeight() ) {
				nextPaths[i] = forgotZero;
			} else {
				nextPaths[i] = forgotOne;
			}
		}
		
		// Output states/paths are input states/paths for the next iteration.
		prevPaths = nextPaths;
	}
	
	cout << "Final Weight     : " << prevPaths[0].getWeight() << endl;

	// Get the best path ending in the state ( 0, 0, 0, 0 ).
	Path selectedPath = prevPaths[0].getPath();
	
	// Copy best path to output (removing trailing bits).
	for ( size_t i = 0; i < selectedPath.size() - R; ++i ) {
		mess_dec.push_back( bitset<K>( selectedPath.at( i ) ) );
	}

	return mess_dec;
}

//////////////////////////////////////////////////////////////////
//                             MAIN                             //
//////////////////////////////////////////////////////////////////

int main()
{
	int NbMot = 12;

	vector< bitset<K> > mess;
	vector< bitset<N> > mess_cod;
	vector< bitset<N> > mess_tra;
	vector< bitset<K> > mess_dec;

	// Random initialization message
	srand( (unsigned)time( NULL ) );
	for ( int i = 0; i < NbMot; ++i )
		mess.push_back( randBitset<K>() );
	for ( int i = 0; i < R; ++i )
		mess.push_back( bitset<K>( 0 ) );

	// Coding of the message => mess_cod
	mess_cod = GSM_code( mess );

	// Simulation of a transmission (errors) => mess_tra
	mess_tra = GSM_transmission( mess_cod );

	// Decoding of the transmitted message => mess_dec
	mess_dec = GSM_decode( mess_tra );

	cout << "Source Message   : ";
	for ( vector<bitset<K> >::iterator it = mess.begin() ; it != mess.end(); ++it )
		cout << ' ' << *it;
	cout << '\n';

	cout << "Coded Message    : ";
	for ( vector<bitset<N> >::iterator it = mess_cod.begin() ; it != mess_cod.end(); ++it )
		cout << ' ' << *it;
	cout << '\n';

	cout << "Received Message : ";
	for ( vector<bitset<N> >::iterator it = mess_tra.begin() ; it != mess_tra.end(); ++it )
		cout << ' ' << *it;
	cout << '\n';

	cout << "Decoded Message  : ";
	for ( vector<bitset<K> >::iterator it = mess_dec.begin() ; it != mess_dec.end(); ++it )
		cout << ' ' << *it;
	cout << '\n';
}

