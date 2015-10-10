#include <stdio.h> // fopen(...), fclose(...), fread(...), fwrite(...), fflush(...), fprintf(...)
#include <string.h> // strncmp(...)
#include <stdlib.h> // exit(...)

#ifdef VS
#	include <malloc.h> // alloca(...)
#endif

// TODO-List:
// ==========
// - Huffman vs Fibonacci
// - Parallel
// - ContainerFormat

using tBool = bool;

using tChar8 = char;

using tNat8  = unsigned char;
using tNat16 = unsigned short;
using tNat32 = unsigned int;
using tNat64 = unsigned long long;

using tInt8  = signed char;
using tInt16 = signed short;
using tInt32 = signed int;
using tInt64 = signed long long;

using tReal32 = float;
using tReal64 = double;

template <typename gRes>
using tFunc0 = gRes();

template <typename gArg, typename gRes>
using tFunc1 = gRes(gArg);

template <typename gArg1, typename gArg2, typename gRes>
using tFunc2 = gRes(gArg1, gArg2);

template <typename gArg1, typename gArg2, typename gArg3, typename gRes>
using tFunc3 = gRes(gArg1, gArg2, gArg3);

using tSize = size_t;
using tStream = FILE;

static tStream* gLogStream = NULL;

#ifndef _countof
#	define _countof(a) (sizeof(a) / sizeof(a[0]))
#endif

#ifdef DEBUG
#	define TEST
#	define TIMER
	
#	define ASSERT(Cond) Debug::Assert((Cond), __FILE__, __LINE__);
	
	namespace Debug {
		
		//=============================================================================
		static inline
		void Assert(
			tBool         Cond,
			const tChar8* FileName,
			tNat32        LineNr
		//=============================================================================
		) {
			if (!Cond) {
				fprintf(gLogStream, "ERROR %s:%d", FileName, LineNr);
				fflush(gLogStream);
				fprintf(gLogStream, "%d", LineNr / 0);
			}
		}
		
	}
#else
#	define ASSERT(Cond)
#endif

//=============================================================================
static inline
tInt32 Abs(
	tInt32 Value
//=============================================================================
) {
	return (Value < 0) ? -Value : Value;
}

//=============================================================================
static inline
tInt32 Sign(
	tInt64 Value
//=============================================================================
) {
	if (Value == 0) { return 0; }
	return (tNat32)(Value >> 63);
}

//=============================================================================
static inline
tInt32 Min(
	tInt32 A,
	tInt32 B
//=============================================================================
) {
	return (A < B) ? A : B;
}

//=============================================================================
static inline
tNat64 Min(
	tNat64 A,
	tNat64 B
//=============================================================================
) {
	return (A < B) ? A : B;
}

//=============================================================================
static inline
tInt32 Max(
	tInt32 A,
	tInt32 B
//=============================================================================
) {
	return (A > B) ? A : B;
}

//=============================================================================
static inline
tNat64 UsedBits(
	tNat64 Nat
//=============================================================================
) {
	auto Bits = 0;
	for (auto Temp = Nat; Temp > 0; Temp >>= 1) {
		Bits += 1;
	}
	return Bits;
}

//=============================================================================
template <typename g>
static inline
void Swap(
	g* A, // IN-OUT
	g* B  // IN OUT
//=============================================================================
) {
	auto Temp = *A;
	*A = *B;
	*B = Temp;
}

template <typename g>
struct tArray sealed {
	tNat64 Size;
	g*     Values;
	
	g& operator[] (const tNat64);
};

//=============================================================================
template <typename g>
inline
g& tArray<g>::operator[] (
	const tNat64 Index
//=============================================================================
) {
	ASSERT(0 <= Index && Index < Size);
	return this->Values[Index];
}

//=============================================================================
template<typename g>
static inline
tArray<g> Take(
	tArray<g> Array,
	tNat32    Count
//=============================================================================
) {
	Array.Size = Min(Array.Size, Count);
	return Array;
}

//=============================================================================
template<typename g>
static inline
tArray<g> Skip(
	tArray<g> Array,
	tNat32    Count
//=============================================================================
) {
	auto Delta = Min(Array.Size, Count);
	Array.Values += Delta;
	Array.Size   -= Delta;
	return Array;
}

//=============================================================================
template<typename g>
static inline
void InPlaceQuickSort(
	tArray<g>            Nodes,
	tFunc2<g, g, tInt32>* Comp
//=============================================================================
) {
	if (Nodes.Size <= 1) {
		return;
	}
	
	auto Pivot  = Nodes[Nodes.Size >> 1];
	auto Low   = (tInt32)0;
	auto Hight = (tInt32)(Nodes.Size - 1);
	
	for (;;) {
		while (Comp(Nodes[Low  ], Pivot) > 0 && Low != Hight) { Low   += 1; }
		while (Comp(Nodes[Hight], Pivot) < 0 && Low != Hight) { Hight -= 1; }
		if (Low == Hight) { break; }
		if (Comp(Nodes[Low], Nodes[Hight]) == 0) {
			Hight -= 1;
		} else {
			Swap(&Nodes[Low], &Nodes[Hight]);
		}
	}
	ASSERT(Low == Hight);
	InPlaceQuickSort(Take(Nodes, Low), Comp);
	InPlaceQuickSort(Skip(Nodes, Low + 1), Comp);
}

#define AsArray(LocalArray) { _countof(LocalArray), LocalArray }

#ifdef TEST
	//=============================================================================
	static inline
	void TestSwap(
	//=============================================================================
	) {
		auto A = (tNat32)5;
		auto B = (tNat32)8;
		
		Swap(&A, &B);
		
		ASSERT(A == 8);
		ASSERT(B == 5);
	}
	
	//=============================================================================
	static inline
	void TestSort(
	//=============================================================================
	) {
		tNat32 Buffer[] = {
			0, // count
			
			1, // count
			1, // in
			1, // out
			
			2, // count
			1, 2, // in
			2, 1, // out
			
			2, // count
			1, 1, // in
			1, 1, // out
			
			2, // count
			2, 1, // in
			2, 1, // out
			
			5, // count
			4, 8, 1, 1, 4, // in
			8, 4, 4, 1, 1, // out
			
			5, // count
			1, 1, 1, 1, 1, // in
			1, 1, 1, 1, 1, // out
			
			9, // count
			4, 3, 8, 2, 4, 2, 8, 23, 1, // in
			23, 8, 8, 4, 4, 3, 2, 2, 1 // out
		};
		tArray<tNat32> Array = AsArray(Buffer);
		auto CountIndex = (tNat32)0;
		auto Comp = (tFunc2<tNat32, tNat32, tInt32>*) [](tNat32 A, tNat32 B) { return (tInt32)A - (tInt32)B; };
		while (CountIndex < Array.Size) {
			auto Count = Array[CountIndex];
			auto T = Take(Skip(Array, CountIndex + 1        ), Count);
			auto R = Take(Skip(Array, CountIndex + 1 + Count), Count);
			InPlaceQuickSort(T, Comp);
			for (auto I = Count; I --> 0;) {
				ASSERT(T[I] == R[I]);
			}
			CountIndex += 2*Count + 1;
		}
	}
#endif

#ifdef TIMER
#	define START_CLOCK(ClockId) Clock::Start(Clock::Id::ClockId);
#	define STOP_CLOCK(ClockId) Clock::Stop(Clock::Id::ClockId);
#	define PRINT_CLOCKS(stream) Clock::print((stream));
	
	namespace Clock {
		
		using tClock = struct {
			tNat64 Begin;
			tNat64 Duration;
			tNat64 Count;
		};
		
		namespace Id {
			enum {
				main,
				DeCode,
				EnCode,
				GetLayerFromBmp16,
				GetLayerFromBmp24,
				GetLayerFromBmp32,
				SplitWithQubPred,
				ComposeWithQubPred,
				ReadLayer,
				WriteLayer,
				HilbertCurve_Next,
				WriteNat_,
				ReadNat_,
				Fib_Write,
				Fib_Read,
				Elias_Write,
				Elias_Read,
				FibToNat,
				NatToFib,
				WriteBit,
				ReadBit,
				ArithmeticBitStream_WriteBit,
				ArithmeticBitStream_ReadBit,
				WriteByte,
				ReadByte,
				WriteByte_fwrite,
				ReadByte_fread,
				PutLayerToBmp32,
				
				_50 = 50,
				_51,
				_52,
				_53,
				_54,
				_55,
				_56,
				_57,
				_58,
				_59,
			};
		}
		
		const tNat32 MaxClocks = 100;
		tClock gClocks[MaxClocks] = {};
		
#		ifdef VS
#			include <intrin.h>
			
			static inline
			//=============================================================================
			tNat64 clock(
			//=============================================================================
			){
				return __rdtsc();
			}
#		else
			static inline
			//=============================================================================
			tNat64 clock(
			//=============================================================================
			) {
				tNat32 lo;
				tNat32 hi;
				__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi) : : "memory");
				return ((tNat64)hi << 32) | lo;
			}
#		endif
		
		//=============================================================================
		static inline
		void Start(
			tNat32 Id
		//=============================================================================
		) {
			gClocks[Id].Begin = clock();
		}
		
		//=============================================================================
		static inline
		void Stop(
			tNat32 Id
		//=============================================================================
		) {
			if (gClocks[Id].Begin) {
				gClocks[Id].Duration += clock() - gClocks[Id].Begin;
				gClocks[Id].Begin = 0;
				gClocks[Id].Count += 1;
			}
		}
		
		//=============================================================================
		static inline
		void print(
			tStream* Stream
		//=============================================================================
		) {
			static char* Names[MaxClocks] = {};
#			define ADD_CLOCK(ClockId) Names[Id::ClockId] = #ClockId;
			ADD_CLOCK(ArithmeticBitStream_ReadBit);
			ADD_CLOCK(ArithmeticBitStream_WriteBit);
			ADD_CLOCK(ComposeWithQubPred);
			ADD_CLOCK(DeCode);
			ADD_CLOCK(EnCode);
			ADD_CLOCK(FibToNat);
			ADD_CLOCK(GetLayerFromBmp16);
			ADD_CLOCK(GetLayerFromBmp24);
			ADD_CLOCK(GetLayerFromBmp32);
			ADD_CLOCK(HilbertCurve_Next);
			ADD_CLOCK(main);
			ADD_CLOCK(NatToFib);
			ADD_CLOCK(PutLayerToBmp32);
			ADD_CLOCK(ReadBit);
			ADD_CLOCK(ReadByte);
			ADD_CLOCK(Fib_Read);
			ADD_CLOCK(Elias_Read);
			ADD_CLOCK(ReadLayer);
			ADD_CLOCK(ReadNat_);
			ADD_CLOCK(SplitWithQubPred);
			ADD_CLOCK(WriteBit);
			ADD_CLOCK(WriteByte);
			ADD_CLOCK(Fib_Write);
			ADD_CLOCK(Elias_Write);
			ADD_CLOCK(WriteLayer);
			ADD_CLOCK(WriteByte_fwrite);
			ADD_CLOCK(WriteNat_);
			ADD_CLOCK(ReadByte_fread);
#			undef ADD_CLOCK
			
			for (auto I = 0; I < MaxClocks; I += 1) {
				auto Count = Clock::gClocks[I].Count;
				if (Count) {
					auto Duration = Clock::gClocks[I].Duration;
					auto Mean = Duration / Count;
					
					auto CountUnit = "";
					if (Count > 10000) { Count /= 1000; CountUnit = "k"; }
					if (Count > 10000) { Count /= 1000; CountUnit = "M"; }
					
					auto DuratuonUnit = "";
					if (Duration > 10000) { Duration /= 1000; DuratuonUnit = "k"; }
					if (Duration > 10000) { Duration /= 1000; DuratuonUnit = "M"; }
					
					auto MeanUnit = "";
					if (Mean > 10000) { Mean /= 1000; MeanUnit = "k"; }
					if (Mean > 10000) { Mean /= 1000; MeanUnit = "M"; }
					
					if (Names[I]) {
						fprintf(Stream, "ClockId %s: ", Names[I]);
					} else {
						fprintf(Stream, "ClockId %d: ", I);
					}
					fprintf(Stream, "%lu", Count);
					fprintf(Stream, "%s ", CountUnit);
					fprintf(Stream, "* %lu", Mean);
					fprintf(Stream, "%s = ", MeanUnit);
					fprintf(Stream, "%lu", Duration);
					fprintf(Stream, "%s\n", DuratuonUnit);
				}
			}
		}
		
	}
#else
#	define START_CLOCK(ClockId)
#	define STOP_CLOCK(ClockId)
#	define PRINT_CLOCKS(stream)
#endif

//=============================================================================
template <typename g>
tArray<g> HeapAlloc(
	tNat64 Count
//=============================================================================
) {
	tArray<g> Result;
	Result.Values = (g*)malloc(Count*sizeof(g));
	Result.Size = Count;
	return Result;
}

//=============================================================================
template <typename g>
void HeapFree(
	tArray<g> Array
//=============================================================================
) {
	free(Array.Values);
	Array.Values = NULL;
	Array.Size = 0;
}

#if 0 // TODO: ByteReaderInterface
using iByteReader = struct {
	void*                 Env;
	tFunc1<void*, tNat32> ReadByte;
	tFunc1<void*, void>   Close;
};

static inline
//=============================================================================
tNat32 ReadByte (
	iByteReader* Reader
//=============================================================================
) {
	return Reader->ReadByte(Reader->Env);
}

static inline
//=============================================================================
void Close (
	iByteReader* Reader
//=============================================================================
) {
	Reader->Close(Reader->Env);
}

using iByteWriter = struct {
	void*                       Env;
	tFunc2<void*, tNat32, void> WriteByte;
	tFunc1<void*, void>         Close;
};

static inline
//=============================================================================
void WriteByte (
	iByteWriter* Writer,
	tNat32       Byte
//=============================================================================
) {
	Writer->WriteByte(Writer->Env, Byte);
}

static inline
//=============================================================================
void Close (
	iByteWriter* Writer
//=============================================================================
) {
	Writer->Close(Writer->Env);
}
#endif

namespace Huffman { // TODO: Huffman
	using tCount_NodeIndex = struct {
		tNat32 Count;
		tNat32 NodeIndex;
	};
	
	// TODO: why i need references ??? (compiler bug?)
	//=============================================================================
	static inline
	void CalculateBitCount(
		tArray<tNat32>& SortedHistogramm,
		tArray<tNat32>& BitCounts         // OUT
	//=============================================================================
	) {
		ASSERT(SortedHistogramm.Size <= (1 << 17))
		
		// init
		auto HistogrammSize = SortedHistogramm.Size;
		auto NodeCount = 2*HistogrammSize - 1;
		
		auto NodeIds     = HeapAlloc<tNat32>(NodeCount);
		auto ParentNodes = HeapAlloc<tNat32>(NodeCount);
		auto Counts      = HeapAlloc<tNat32>(NodeCount);
		
		for (auto NodeId = NodeCount; NodeId --> 0;) {
			NodeIds[NodeId] = (tNat32)NodeId;
			ParentNodes[NodeId] = 0;
			Counts[NodeId] = (NodeId < HistogrammSize) ? SortedHistogramm[NodeId] : 0;
		}
		
		{ // build tree
			auto I = (tNat32)0;
			auto NewNodeId = (tNat32)HistogrammSize;
			while (I <= NodeCount - 2) {
				ASSERT(NewNodeId <= NodeCount)
				auto NodeIdA = NodeIds[I + 0];
				auto NodeIdB = NodeIds[I + 1];
				ParentNodes[NodeIdA] = NewNodeId;
				ParentNodes[NodeIdB] = NewNodeId;
				
				Counts[NewNodeId] = Counts[NodeIdA] + Counts[NodeIdB];
				
				// insert sort, one round
				for (auto J = NewNodeId; J > I + 2; J -= 1) {
					if (Counts[NodeIds[J + 0]] > Counts[NodeIds[J - 1]]) { break; }
					Swap(&NodeIds[J + 0], &NodeIds[J - 1]);
				}
				
				I += 2;
				NewNodeId += 1;
			}
		}
		
		// get bit count
		Counts[NodeCount - 1] = 0; // root
		for (auto I = NodeCount - 1; I --> 0;) {
			auto NodeId       = NodeIds[I];
			auto ParentNodeId = ParentNodes[NodeId];
			Counts[NodeId]    = Counts[ParentNodeId] + 1;
			if (NodeId < HistogrammSize) {
				BitCounts[NodeId] = Counts[NodeId];
			}
		}
		
		// clean up
		HeapFree(ParentNodes);
		HeapFree(NodeIds);
		HeapFree(Counts);
	}
	

	//=============================================================================
	static inline
	void CalculateBits(
		tArray<tNat32>& BitCounts,
		tArray<tNat32>& BitsArray  // OUT
	//=============================================================================
	) {
		auto Bits = 0;
		auto BitCount = 0;
		for (auto NodeId = BitCounts.Size; NodeId --> 0;) {
			auto NewBitCount = BitCounts[NodeId];
			Bits <<= NewBitCount - BitCount;
			BitCount = NewBitCount;
			
			BitsArray[NodeId] = Bits;
			
			Bits += 1;
		}
	}
	
#	ifdef TEST
		//=============================================================================
		static inline
		void TestCalcBitCount(
		//=============================================================================
		) {
			tNat32 Histogram[]    = { 1, 1, 4, 4, 8, 20, 20, 20 };
			//                         \2/  |  |  |   |   |   |
			//                           \6/   |  |   |   |   |
			//                             \10/   |   |   |   |
			//                                \18/    |   |   |
			//                                   \ 38/    |   |
			//                                       \    \40/
			//                                        \78/
			tNat32 BitCountsRef[] = { 6, 6, 5, 4, 3,  2,  2,  2 };
			tNat32 BitCountsRes[] = { 0, 0, 0, 0, 0,  0,  0,  0 };
			tNat32 BitsRes[] = {        0,        0,       0,      0,     0,    0,    0,    0 };
			tNat32 BitsRef[] = { 0b111111, 0b111110, 0b11110, 0b1110, 0b110, 0b10, 0b01, 0b00 };
			
			tArray<tNat32> InHistogramm = AsArray(Histogram);
			tArray<tNat32> OutCounts    = AsArray(BitCountsRes);
			tArray<tNat32> OutBits      = AsArray(BitsRes);
			
			CalculateBitCount(InHistogramm, OutCounts);
			for (auto I = _countof(Histogram); I --> 0;) {
				ASSERT(BitCountsRes[I] == BitCountsRef[I]);
			}
			
			CalculateBits(OutCounts, OutBits);
			for (auto I = _countof(Histogram); I --> 0;) {
				ASSERT(BitsRes[I] == BitsRef[I]);
			}
		}
		
		//=============================================================================
		static inline
		void Test(
		//=============================================================================
		) {
			TestSwap();
			TestSort();
			TestCalcBitCount();
		}
#	endif
}

namespace BufferdStream {
	using tReader = struct {
		tNat8    Buffer[512];
		tStream* Stream;
		tNat32   Pos;
		tNat32   _; // aligment
	};
	
	using tWriter = struct {
		tNat8    Buffer[512];
		tStream* Stream;
		tNat32   Pos;
		tNat32   _; // aligment
	};
	
	static inline
	//=============================================================================
	void Init(
		tReader* BufferdStream,
		tStream* Stream
	//=============================================================================
	) {
		BufferdStream->Pos = 0;
		BufferdStream->Stream = Stream;
	}
	
	static inline
	//=============================================================================
	void Init(
		tWriter* BufferdStream,
		tStream* Stream
	//=============================================================================
	) {
		BufferdStream->Pos = 0;
		BufferdStream->Stream = Stream;
	}
	
	static inline
	//=============================================================================
	void WriteByte(
		tWriter* Stream,
		tNat8    Byte
	//=============================================================================
	) {
//		START_CLOCK(WriteByte);
		
		const tNat32 BufferSize = sizeof(Stream->Buffer);
		auto Pos = Stream->Pos;
		Stream->Buffer[Pos] = Byte;
		Pos = (Pos + 1) % BufferSize;
		if (Pos == 0) {
//			START_CLOCK(WriteByte_fwrite);
			fwrite(Stream->Buffer, BufferSize, 1, Stream->Stream);
//			STOP_CLOCK(WriteByte_fwrite);
		}
		Stream->Pos = Pos;
		
//		STOP_CLOCK(WriteByte);
	}
	
	static inline
	//=============================================================================
	tInt32 ReadByte(
		tReader* Stream
	//=============================================================================
	) {
//		START_CLOCK(ReadByte);
		
		const tNat32 BufferSize = sizeof(Stream->Buffer);
		auto Pos = Stream->Pos;
		if (Pos == 0) {
//			START_CLOCK(ReadByte_fread);
			auto Bytes = (tNat32)fread(Stream->Buffer, 1, BufferSize, Stream->Stream);
//			STOP_CLOCK(ReadByte_fread);
			if (Bytes == 0) {
				STOP_CLOCK(ReadByte);
				return -1; // TODO: Fehlerbehandlung
			}
			if (Bytes != BufferSize) {
				auto Des = BufferSize;
				auto Src = Bytes;
				while (Src > 0) {
					Des -= 1;
					Src -= 1;
					Stream->Buffer[Des] = Stream->Buffer[Src];
				}
				Pos = BufferSize - Bytes;
			}
		}
		Stream->Pos = (Pos + 1) % BufferSize;
		auto Result =  Stream->Buffer[Pos];
		
//		STOP_CLOCK(ReadByte);
		return Result;
	}
	
	static inline
	//=============================================================================
	void Close(
		tReader* Stream
	//=============================================================================
	) {
		Stream->Pos = 0;
		Stream->Stream = NULL;
	}
	
	static inline
	//=============================================================================
	void Close(
		tWriter* Stream
	//=============================================================================
	) {
		const tNat32 BufferSize = sizeof(Stream->Buffer);
//		START_CLOCK(WriteByte_fwrite);
		fwrite(Stream->Buffer, Stream->Pos, 1, Stream->Stream);
//		STOP_CLOCK(WriteByte_fwrite);
		Stream->Pos = 0;
		Stream->Stream = NULL;
	}
}

using iBitReader = struct {
	tFunc1<void*, tNat32>* ReadBit;
	tFunc1<void*, void>*   Close;
	void*                  Context;
};

//=============================================================================
static inline
tNat32 ReadBit(
	iBitReader* BitReader
//=============================================================================
) {
	return BitReader->ReadBit(BitReader->Context);
}

//=============================================================================
static inline
void Close(
	iBitReader* BitReader
//=============================================================================
) {
	BitReader->Close(BitReader->Context);
	*BitReader = {};
}

using iBitWriter = struct {
	tFunc2<void*, tNat32, void>* WriteBit;
	tFunc1<void*, void>*         Close;
	void*                        Context;
};

//=============================================================================
static inline
void WriteBit(
	iBitWriter* BitWriter,
	tNat32      Bit
//=============================================================================
) {
	BitWriter->WriteBit(BitWriter->Context, Bit);
}

//=============================================================================
static inline
void Close(
	iBitWriter* BitWriter
//=============================================================================
) {
	BitWriter->Close(BitWriter->Context);
	*BitWriter = {};
}

namespace BitStream {
	using tReader = struct {
		BufferdStream::tReader* Stream;
		tNat32                  Bits;
		tNat32                  Byte;
	};
	
	using tWriter = struct {
		BufferdStream::tWriter* Stream;
		tNat32                  Bits;
		tNat32                  Byte;
	};
	
	//=============================================================================
	static inline
	void Init(
		tWriter*                BitStream,
		BufferdStream::tWriter* Stream
	//=============================================================================
	) {
		*BitStream = {};
		BitStream->Stream = Stream;
	}
	
	//=============================================================================
	static inline
	void Close(
		tWriter* BitStream
	//=============================================================================
	) {
		if (BitStream->Bits) {
			BufferdStream::WriteByte(BitStream->Stream, (tNat8)BitStream->Byte);
		}
		*BitStream = {};
	}
	
	//=============================================================================
	static inline
	void WriteBit(
		tWriter* BitStream,
		tNat32   Bit
	//=============================================================================
	) {
//		START_CLOCK(WriteBit);
		
		ASSERT((Bit & ~1) == 0);
		
		auto Bits = BitStream->Bits;
		
		if (Bits == 0) {
			BitStream->Byte = Bit;
		} else {
			for (auto I = Bits; I --> 0;) { Bit <<= 1; }
			BitStream->Byte |= Bit;
		}
		
		Bits = (Bits + 1) & 7;
		if (Bits == 0) {
			BufferdStream::WriteByte(BitStream->Stream, (tNat8)BitStream->Byte);
		}
		BitStream->Bits = Bits;
		
//		STOP_CLOCK(WriteBit);
	}
	
	//=============================================================================
	static inline
	void Init(
		tReader*                BitStream,
		BufferdStream::tReader* Stream
	//=============================================================================
	) {
		*BitStream = {};
		BitStream->Stream = Stream;
	}
	
	//=============================================================================
	static inline
	void Close(
		tReader* BitStream
	//=============================================================================
	) {
		*BitStream = {};
	}
	
	//=============================================================================
	static inline
	tNat32 ReadBit(
		tReader* BitStream
	//=============================================================================
	) {
//		START_CLOCK(ReadBit);
		
		auto Bits = BitStream->Bits;
		auto Byte = 0;
		if (Bits == 0) {
			Byte = BufferdStream::ReadByte(BitStream->Stream);
		} else {
			Byte = BitStream->Byte >> 1;
		}
		
		Bits = (Bits + 1) & 7;
		
		BitStream->Bits = Bits;
		BitStream->Byte = Byte;
		
//		STOP_CLOCK(ReadBit);
		return Byte & 1;
	}
	
	//=============================================================================
	static inline
	iBitWriter GetInterface(
		tWriter* BitStream
	//=============================================================================
	) {
		iBitWriter Interface = {};
		Interface.Context  = BitStream;
		Interface.WriteBit = (tFunc2<void*, tNat32, void>*)WriteBit;
		Interface.Close    = (tFunc1<void*, void>*)(tFunc1<tWriter*, void>*)Close;
		return Interface;
	}
	
	//=============================================================================
	static inline
	iBitReader GetInterface(
		tReader* BitStream
	//=============================================================================
	) {
		iBitReader Interface = {};
		Interface.Context = BitStream;
		Interface.ReadBit = (tFunc1<void*, tNat32>*)ReadBit;
		Interface.Close   = (tFunc1<void*, void>*)(tFunc1<tReader*, void>*)Close;
		return Interface;
	}
}

namespace ArithmeticBitStream {
	using tWriter = struct {
		iBitWriter* BitStream;
		tNat32 UnknowBits;
		tNat32 Min;
		tNat32 Range;
		tNat32 Count0;
		tNat32 Count;
		tNat32 Count0_Old;
	};
	
	using tReader = struct {
		iBitReader* BitStream;
		tNat32 UnknowBits;
		tNat32 Value;
		tNat32 Min;
		tNat32 Range;
		tNat32 Count0;
		tNat32 Count;
		tNat32 Count0_Old;
		tInt32 _; // 64 bit alignment
	};
	
	//=============================================================================
	static inline
	void Init(
		tWriter*    Writer,
		iBitWriter* BitStream
	//=============================================================================
	) {
		*Writer = {};
		Writer->BitStream  = BitStream;
		Writer->Range      = 1 << 16;
		Writer->Count0     = 1 << 7;
		Writer->Count0_Old = 1 << 7;
	}
	
	//=============================================================================
	static inline
	void Init(
		tReader*    Reader,
		iBitReader* BitStream
	//=============================================================================
	) {
		*Reader = {};
		Reader->BitStream  = BitStream;
		Reader->Range      = 1 << 16;
		Reader->Count0     = 1 << 7;
		Reader->Count0_Old = 1 << 7;
	}
	
	//=============================================================================
	static inline
	void WriteBit(
		tWriter* Stream,
		tNat32   Bit
	//=============================================================================
	) {
//		START_CLOCK(ArithmeticBitStream_WriteBit);
		
		const tNat32 Full = 1 << 16;
		const tNat32 Half = 1 << 15;
		const tNat32 Quad = 1 << 14;
		
		auto Delta = (Stream->Count0_Old * Stream->Range) >> 8;
		
		if (Bit == 0) {
			Stream->Range  = Delta;
		} else {
			Stream->Min   += Delta;
			Stream->Range -= Delta;
		}
		
		for(;;) {
			if ((Stream->Min + Stream->Range) - 1 < Half) {
				Stream->Range <<= 1;
				Stream->Min    -= 0;
				Stream->Min   <<= 1;
				
				WriteBit(Stream->BitStream, 0);
				while (Stream->UnknowBits) {
					Stream->UnknowBits -= 1;
					WriteBit(Stream->BitStream, 1);
				}
			} else if (Stream->Min >= Half) {
				Stream->Range <<= 1;
				Stream->Min    -= Half;
				Stream->Min   <<= 1;
				
				WriteBit(Stream->BitStream, 1);
				while (Stream->UnknowBits) {
					Stream->UnknowBits -= 1;
					WriteBit(Stream->BitStream, 0);
				}
			} else if (
				Stream->Min >= Quad &&
				(Stream->Min + Stream->Range) - 1 < 3*Quad
			) {
				Stream->Range <<= 1;
				Stream->Min    -= Quad;
				Stream->Min   <<= 1;
				
				Stream->UnknowBits += 1;
			} else {
				break;
			}
		}
		
		if (Bit == 0) {
			Stream->Count0 += 1;
		}
		Stream->Count = (Stream->Count + 1) & 0xFF;
		if (Stream->Count == 0) {
			Stream->Count0   >>= 1;
			Stream->Count0_Old = (tNat8)Stream->Count0;
		}
		
//		START_CLOCK(ArithmeticBitStream_WriteBit);
	}
	
	//=============================================================================
	static inline
	tNat32 ReadBit(
		tReader* Stream
	//=============================================================================
	) {
//		START_CLOCK(ArithmeticBitStream_ReadBit);
		
		const tNat32 Full = 1 << 16;
		const tNat32 Half = 1 << 15;
		const tNat32 Quad = 1 << 14;
		
		auto Bit = (tNat32)0;
		
		if (Stream->UnknowBits) {
			Bit = Stream->Value;
			Stream->UnknowBits -= 1;
		} else {
			for (;;) {
				if ((Stream->Min + Stream->Range) - 1 < Half) {
					Bit = 0;
					Stream->Value = 1;
					break;
				} else if (Stream->Min >= Half) {
					Bit = 1;
					Stream->Value = 0;
					break;
				} else if (
					Stream->Min >= Quad &&
					(Stream->Min + Stream->Range) - 1 < 3*Quad
				) {
					Stream->UnknowBits += 1;
					Stream->Range <<= 1;
					Stream->Min    -= Quad;
					Stream->Min   <<= 1;
				} else {
					if (ReadBit(Stream->BitStream) == 0) {
						Stream->Range <<= 1;
						Stream->Min    -= 0;
						Stream->Min   <<= 1;
					} else {
						Stream->Range <<= 1;
						Stream->Min    -= Half;
						Stream->Min   <<= 1;
					}
				}
			}
		}
		
		if (Bit == 0) {
			Stream-> Count0 += 1;
		}
		Stream->Count = (Stream->Count + 1) & 0xFF;
		if (Stream->Count == 0) {
			Stream->Count0_Old = (tNat8)Stream->Count0;
			Stream->Count0   >>= 1;
		}
		
		auto Delta = (Stream->Count0_Old * Stream->Range) >> 8;
		
		if (Bit == 0) {
			Stream->Range  = Delta;
		} else {
			Stream->Min   += Delta;
			Stream->Range -= Delta;
		}
		
//		STOP_CLOCK(ArithmeticBitStream_ReadBit);
		return Bit;
	}
	
	//=============================================================================
	static inline
	void Close(
		tWriter* Stream
	//=============================================================================
	) {
		Close(Stream->BitStream);
		*Stream = {};
	}
	
	//=============================================================================
	static inline
	void Close(
		tReader* Stream
	//=============================================================================
	) {
		Close(Stream->BitStream);
		*Stream = {};
	}
	
	//=============================================================================
	static inline
	iBitReader GetInterface(
		tReader* Reader
	//=============================================================================
	) {
		iBitReader Interface = {};
		Interface.Context = (void*)Reader;
		Interface.ReadBit = (tFunc1<void*, tNat32>*)(tFunc1<tReader*, tNat32>*)ReadBit;
		Interface.Close   = (tFunc1<void*, void>*)(tFunc1<tReader*, void>*)Close;
		return Interface;
	}
	
	//=============================================================================
	static inline
	iBitWriter GetInterface(
		tWriter* Writer
	//=============================================================================
	) {
		iBitWriter Interface = {};
		Interface.Context  = (void*)Writer;
		Interface.WriteBit = (tFunc2<void*, tNat32, void>*)(tFunc2<tWriter*, tNat32, void>*)WriteBit;
		Interface.Close    = (tFunc1<void*, void>*)(tFunc1<tWriter*, void>*)Close;
		return Interface;
	}
}

//=============================================================================
static inline
void WriteNat_(
	iBitWriter* BitWriter,
	tNat64      Value
//=============================================================================
) {
//	START_CLOCK(WriteNat_);
	
	for (auto I = (tNat32)16; I --> 0;) {
		auto Bit = (Value >> I) & 1;
		WriteBit(BitWriter, (tNat32)Bit);
	}
	
//	STOP_CLOCK(WriteNat_);
}

//=============================================================================
static inline
tNat64 ReadNat_(
	iBitReader* BitReader
//=============================================================================
) {
//	START_CLOCK(ReadNat_);
	
	auto Result = (tNat64)0;
	for (auto I = 16; I --> 0;) {
		Result <<= 1;
		Result |= ReadBit(BitReader);
	}
	
//	STOP_CLOCK(ReadNat_);
	return Result;
}

namespace EliasGammaCode {
	
	//=============================================================================
	static inline
	void Write(
		iBitWriter* BitWriter,
		tNat64      Value
	//=============================================================================
	) {
//		START_CLOCK(Elias_Write);
		
		Value += 1;
		
		auto DataBits = 1;
		for (auto Temp = Value >> 1; Temp > 0; Temp >>= 1) {
			WriteBit(BitWriter, 0);
			DataBits += 1;
		}
		
		while (DataBits > 0) {
			DataBits -= 1;
			WriteBit(BitWriter, (Value >> DataBits) & 1);
		}
		
//		STOP_CLOCK(Elias_Write);
	}
	
	//=============================================================================
	static inline
	tNat64 Read(
		iBitReader* BitReader
	//=============================================================================
	) {
//		START_CLOCK(Elias_Read);
		
		auto DataBits = 1;
		auto Bit = ReadBit(BitReader);
		while (Bit == 0) {
			DataBits += 1;
			Bit = ReadBit(BitReader);
		}
		
		auto Value = 1;
		while (DataBits -= 1) {
			Value <<= 1;
			Value  |= ReadBit(BitReader);
		}
		
//		STOP_CLOCK(Elias_Read);
		return Value - 1;
	}
	
}

namespace FibonacciCode {
	static tNat64 gFibArray[50];
	
	//=============================================================================
	static inline
	void Init(
		tNat64* FibArray,
		tNat32  MaxFibBits
	//=============================================================================
	) {
		auto F1 = (tNat64)1;
		auto F2 = (tNat64)1;

		ASSERT(MaxFibBits == 50)
		
		for (auto FibBits = (tNat32)0; FibBits <= MaxFibBits; FibBits += 1) {
			FibArray[FibBits] = F2;
			
			auto F = F1 + F2;
			
			F1 = F2;
			F2 = F;
		}
	}

	//=============================================================================
	static inline
	tNat64 NatToFib(
		tNat64* FibArray,
		tNat64  Nat
	//=============================================================================
	) {
//		START_CLOCK(NatToFib);
		
		auto FibBits = UsedBits(Nat);
		FibBits += FibBits >> 1;
		
		auto Fib = (tNat64)0;
		for (; FibBits > 0; FibBits -= 1) {
			auto FibI = FibArray[FibBits-1];
			auto Bit = Nat >= FibI;
			if (Bit) {
				Nat -= FibI;
			}
			Fib <<= 1;
			Fib |= (tNat64)Bit;
		}
		
//		STOP_CLOCK(NatToFib);
		return Fib;
	}
	
	//=============================================================================
	static inline
	tNat64 FibToNat(
		tNat64 Fib
	//=============================================================================
	) {
//		START_CLOCK(FibToNat);
		
		auto Nat = (tNat64)0;
		
		auto F1 = (tNat64)0;
		auto F2 = (tNat64)1;
		for (auto I = 0; Fib; I += 1, Fib >>= 1) {
			auto F = F1 + F2;
			
			F1 = F2;
			F2 = F;
			
			if (Fib & 1) {
				Nat += F;
			}
		}
		
//		STOP_CLOCK(FibToNat);
		return Nat;
	}
	
	//=============================================================================
	static inline
	void Write(
		iBitWriter* BitWriter,
		tNat64      Value
	//=============================================================================
	) {
//		START_CLOCK(Fib_Write);
		auto Fib = NatToFib(gFibArray, Value + 1);
		
		while (Fib) {
			WriteBit(BitWriter, Fib & 1);
			Fib >>= 1;
		}
		
		WriteBit(BitWriter, 1);
		
//		STOP_CLOCK(Fib_Write);
	}
	
	//=============================================================================
	static inline
	tNat64 Read(
		iBitReader* BitReader
	//=============================================================================
	) {
//		START_CLOCK(Fib_Read);
		
		auto Fib = (tNat64)ReadBit(BitReader);
		for (auto Pos = 1;; Pos += 1) {
			auto Bit = (tNat64)ReadBit(BitReader);
			if (Bit & (Fib >> (Pos - 1))) { break; }
			Fib |= Bit << Pos;
		}
		
		auto Result = FibToNat(Fib) - 1;
		
//		STOP_CLOCK(Fib_Read);
		return Result;
	}
	
#	ifdef TEST
		//=============================================================================
		static inline
		void Test(
		//=============================================================================
		) {
			Init(gFibArray, _countof(gFibArray));
			for (auto Nat = (tNat32)1; Nat < 10000; Nat += 1) {
				auto Fib = NatToFib(gFibArray, Nat);
				ASSERT(FibToNat(Fib) == Nat);
			}
		}
#	endif
}

namespace HaarWavelet {
	using tQuad = struct {
		tInt16 A;
		tInt16 B;
		tInt16 C;
		tInt16 D;
	};
	
	//=============================================================================
	static inline
	tQuad EnCode(
		tQuad Src
	//=============================================================================
	) {
		tQuad Res;
#		if 1
			Res.A = (Src.A + Src.B + Src.C + Src.D + 2)>>2;
			Res.B =  Src.A - Src.B + Src.C - Src.D;
			Res.C =  Src.A + Src.B - Src.C - Src.D;
			Res.D =  Src.A - Src.B - Src.C + Src.D;
#		else
			Res.A = (Src.A + Src.B + Src.C + Src.D + 2)>>2;
			Res.B =  (Src.A - Src.B + Src.C - Src.D) >> 1;
			Res.C =  (Src.A + Src.B - Src.C - Src.D) >> 1;
			Res.D =  (Src.A - Src.B - Src.C + Src.D) >> 1;
#		endif
		return Res;
	}
	
	//=============================================================================
	static inline
	tQuad DeCode(
		tQuad Src
	//=============================================================================
	) {
		tQuad Res;
#		if 1
			Res.A = ((Src.A<<2) + 1 + Src.B + Src.C + Src.D)>>2;
			Res.B = ((Src.A<<2) + 1 - Src.B + Src.C - Src.D)>>2;
			Res.C = ((Src.A<<2) + 1 + Src.B - Src.C - Src.D)>>2;
			Res.D = ((Src.A<<2) + 1 - Src.B - Src.C + Src.D)>>2;
#		else
			Res.A = ((Src.A<<1) + 1 + Src.B + Src.C + Src.D)>>1;
			Res.B = ((Src.A<<1) + 1 - Src.B + Src.C - Src.D)>>1;
			Res.C = ((Src.A<<1) + 1 + Src.B - Src.C - Src.D)>>1;
			Res.D = ((Src.A<<1) + 1 - Src.B - Src.C + Src.D)>>1;
#		endif
		return Res;
	}
	
#	ifdef TEST
		//=============================================================================
		static inline
		void Test(
		//=============================================================================
		) {
			tQuad Src;
			for (Src.A = -4; Src.A < 4; Src.A += 1) {
				for (Src.B = -4; Src.B < 4; Src.B += 1) {
					for (Src.C = -4; Src.C < 4; Src.C += 1) {
						for (Src.D = -4; Src.D < 4; Src.D += 1) {
							auto Temp = HaarWavelet::EnCode(Src);
							auto Res = HaarWavelet::DeCode(Temp);
							ASSERT(Res.A == Src.A);
							ASSERT(Res.B == Src.B);
							ASSERT(Res.C == Src.C);
							ASSERT(Res.D == Src.D);
						}
					}
				}
			}
		}
#	endif
}

namespace Layer {
	using tLevel = struct {
		tNat32 SizeXLeft;
		tNat32 SizeXRight;       //    Left  Right
		tNat32 SizeYTop;         //   +-----+---+
		tNat32 SizeYBottom;      //   |  S  | H | Top
		tInt16* S; // Summe      //   |     |   |
		tInt16* H; // Horizontal //   +-----+---+
		tInt16* V; // Vertical   //   |  V  | D | Bottom
		tInt16* D; // Diagonal   //   +-----+---+
		tNat32 Pitch;            // Pitch >= SizeXLeft >= SizeXRight; SizeYTop >= SizeYBottom
		tInt32 _; // 64 bit alignment
	};
	
	//=============================================================================
	static
	tNat32 Init(
		tLevel* Layers,
		tNat32  SizeX,
		tNat32  SizeY,
		tInt16* Buffer,
		tInt32  BufferSize
	//=============================================================================
	) {
		auto InitLevel = 0u;
		{ // Berechne InitLevel + Initalisierung vom InitLevel
			auto X = SizeX;
			auto Y = SizeY;
			while (X > 1 || Y > 1) {
				X = (tNat16)((X + 1) >> 1);
				Y = (tNat16)((Y + 1) >> 1);
				InitLevel += 1;
			}
			auto Layer = &Layers[InitLevel];
			
			Layer->Pitch = (tNat32)((SizeX + 7) & ~7);
			Layer->SizeXLeft = SizeX;
			Layer->SizeXRight = 0;
			Layer->SizeYTop = SizeY;
			Layer->SizeYBottom = 0;
			
			Layer->S = Buffer;
			Buffer += Layer->Pitch * Layer->SizeYTop;
			BufferSize -= Layer->Pitch * Layer->SizeYTop;
			Layer->H = NULL;
			Layer->V = NULL;
			Layer->D = NULL;
		}
		
		ASSERT(BufferSize > 0)
		
		{ // Initalisierung restlicher Levels
			auto X = (tNat32)SizeX;
			auto Y = (tNat32)SizeY;
			for (auto Level = (tInt32)InitLevel - 1; Level >= 0; Level -= 1) {
				auto Layer = &Layers[Level];
				
				Layer->SizeXRight  = (X + 0) >> 1;
				Layer->SizeYBottom = (Y + 0) >> 1;
				Layer->SizeXLeft   = (X + 1) >> 1;
				Layer->SizeYTop    = (Y + 1) >> 1;
				
				Layer->Pitch = (Layer->SizeXLeft + 7) & ~7;
				
				Layer->S    = Buffer;
				Buffer     += Layer->Pitch * Layer->SizeYTop;
				BufferSize -= Layer->Pitch * Layer->SizeYTop;
				ASSERT(BufferSize > 0)
				
				Layer->H    = Buffer;
				Buffer     += Layer->Pitch * Layer->SizeYTop;
				BufferSize -= Layer->Pitch * Layer->SizeYTop;
				ASSERT(BufferSize > 0)
				
				Layer->V    = Buffer;
				Buffer     += Layer->Pitch * Layer->SizeYBottom;
				BufferSize -= Layer->Pitch * Layer->SizeYBottom;
				ASSERT(BufferSize > 0)
				
				Layer->D    = Buffer;
				Buffer     += Layer->Pitch * Layer->SizeYBottom;
				BufferSize -= Layer->Pitch * Layer->SizeYBottom;
				ASSERT(BufferSize > 0)
				
				X = Layer->SizeXLeft;
				Y = Layer->SizeYTop;
			}
		}
		return InitLevel;
	}
	
	//=============================================================================
	static inline
	tInt32 DeltaQubicWavelet(
		tInt32 Last2,
		tInt32 Last1,
		tInt32 Curr,
		tInt32 Next1,
		tInt32 Next2
	//=============================================================================
	) {
		//auto c = 7; // beste Ganzzahl
		//return ((c+4)*(Last1 - Next1) - 2*(Last2 - Next2) /*+ c*/) / (2*c);
		
		return (6*(Last1 - Next1) - Last2 + Next2 + 4) >> 3; // entspricht c = 8; beste Zweierpotenz
		//return (Last1 - Next1 + 1) >> 1; // entspricht c -> INF; linear
		//return 0; // OFF
	}
	
	//=============================================================================
	static inline
	tInt32 DeltaQubicWaveletPre2(
		tInt32 Curr,
		tInt32 Next1,
		tInt32 Next2
	//=============================================================================
	) {
		auto Last2 = 2*Curr - Next2;
		auto Last1 = 2*Curr - Next1;
		return DeltaQubicWavelet(Last2, Last1, Curr, Next1, Next2);
	}
	
	//=============================================================================
	static inline
	tInt32 DeltaQubicWaveletPre1(
		tInt32 Last1,
		tInt32 Curr,
		tInt32 Next1,
		tInt32 Next2
	//=============================================================================
	) {
		auto Last2 = 2*Last1 - Curr;
		return DeltaQubicWavelet(Last2, Last1, Curr, Next1, Next2);
	}
	
	//=============================================================================
	static inline
	tInt32 DeltaQubicWaveletPost1(
		tInt32 Last2,
		tInt32 Last1,
		tInt32 Curr,
		tInt32 Next1
	//=============================================================================
	) {
		auto Next2 = 2*Next1 - Curr;
		return DeltaQubicWavelet(Last2, Last1, Curr, Next1, Next2);
	}
	
	//=============================================================================
	static inline
	tInt32 DeltaQubicWaveletPost2(
		tInt32 Last2,
		tInt32 Last1,
		tInt32 Curr
	//=============================================================================
	) {
		auto Next2 = 2*Curr - Last2;
		auto Next1 = 2*Curr - Last1;
		return DeltaQubicWavelet(Last2, Last1, Curr, Next1, Next2);
	}
	
	const tInt32 PredDist = 2;
	const tInt32 DFactor = 2;// shift 2
	
	//=============================================================================
	static
	void SplitWithQubPred(
		tLevel* SrcLevel,
		tLevel* DesLevel
	//=============================================================================
	) {
		START_CLOCK(SplitWithQubPred);
		
		auto SizeX    = (tInt32)SrcLevel->SizeXLeft;
		auto SizeY    = (tInt32)SrcLevel->SizeYTop;
		auto SrcPitch = (tInt32)SrcLevel->Pitch;
		auto DesPitch = (tInt32)DesLevel->Pitch;
		auto Src      = SrcLevel->S;
		auto DesS     = DesLevel->S;
		auto DesH     = DesLevel->H;
		auto DesV     = DesLevel->V;
		auto DesD     = DesLevel->D;
		
		auto SizeX0 = (tInt32)(SizeX + 1) >> 1;
		auto SizeX1 = (tInt32)(SizeX + 0) >> 1;
		auto SizeY0 = (tInt32)(SizeY + 1) >> 1;
		auto SizeY1 = (tInt32)(SizeY + 0) >> 1;
		
		auto Last2_S = 0;
		auto Last2_V = 0;
		auto Last1_S = 0;
		auto Last1_V = 0;
		auto Curr0_S = 0;
		auto Curr0_V = 0;
		auto Next1_S = 0;
		auto Next1_V = 0;
		auto Next2_S = 0;
		auto Next2_V = 0;
		
		for (auto Y = -PredDist; Y < SizeY0; Y++) {
			auto Y2 = Y << 1;
			for (auto X = -PredDist; X < SizeX0; X++) {
				Last2_S = Last1_S;
				Last2_V = Last1_V;
				Last1_S = Curr0_S;
				Last1_V = Curr0_V;
				Curr0_S = Next1_S;
				Curr0_V = Next1_V;
				Next1_S = Next2_S;
				Next1_V = Next2_V;
				if (X + PredDist < SizeX0) {
					if (Y + PredDist < SizeY0) {
						auto X2 = X << 1;
						HaarWavelet::tQuad Arg;
						HaarWavelet::tQuad Res;
						if (Y2 + 2*PredDist == SizeY - 1) {
							if (X2 + 2*PredDist  == SizeX - 1) {
								Arg.A = Src[(tNat32)((X2 + 2*PredDist + 0) + (Y2 + 2*PredDist + 0) * SrcPitch)];
								Arg.B = Arg.A;
								Arg.C = Arg.A;
								Arg.D = Arg.A;
								
								Res = HaarWavelet::EnCode(Arg);
								
								DesS[(X + PredDist) + (Y + PredDist) * DesPitch] = Res.A;
							} else {
								Arg.A = Src[(tNat32)((X2 + 2*PredDist + 0) + (Y2 + 2*PredDist + 0) * SrcPitch)];
								Arg.B = Src[(tNat32)((X2 + 2*PredDist + 1) + (Y2 + 2*PredDist + 0) * SrcPitch)];
								Arg.C = Arg.A;
								Arg.D = Arg.B;
								
								Res = HaarWavelet::EnCode(Arg);
								
								DesS[(X + PredDist) + (Y + PredDist) * DesPitch] = Res.A;
								DesH[(X + PredDist) + (Y + PredDist) * DesPitch] = Res.B;
							}
						} else {
							if (X2 + 2*PredDist == SizeX - 1) {
								Arg.A = Src[(tNat32)((X2 + 2*PredDist + 0) + (Y2 + 2*PredDist + 0) * SrcPitch)];
								Arg.B = Arg.A;
								Arg.C = Src[(tNat32)((X2 + 2*PredDist + 0) + (Y2 + 2*PredDist + 1) * SrcPitch)];
								Arg.D = Arg.C;
								
								Res = HaarWavelet::EnCode(Arg);
								
								DesS[(X + PredDist) + (Y + PredDist) * DesPitch] = Res.A;
								DesV[(X + PredDist) + (Y + PredDist) * DesPitch] = Res.C;
							} else {
								Arg.A = Src[(tNat32)((X2 + 2*PredDist + 0) + (Y2 + 2*PredDist + 0) * SrcPitch)];
								Arg.B = Src[(tNat32)((X2 + 2*PredDist + 1) + (Y2 + 2*PredDist + 0) * SrcPitch)];
								Arg.C = Src[(tNat32)((X2 + 2*PredDist + 0) + (Y2 + 2*PredDist + 1) * SrcPitch)];
								Arg.D = Src[(tNat32)((X2 + 2*PredDist + 1) + (Y2 + 2*PredDist + 1) * SrcPitch)];
								
								Res = HaarWavelet::EnCode(Arg);
								
								DesS[(X + PredDist) + (Y + PredDist) * DesPitch] = Res.A;
								DesH[(X + PredDist) + (Y + PredDist) * DesPitch] = Res.B;
								DesV[(X + PredDist) + (Y + PredDist) * DesPitch] = Res.C;
								DesD[(X + PredDist) + (Y + PredDist) * DesPitch] = Res.D;
							}
						}
					}
					
					if (Y < 0 || (SizeX & SizeY) < 8) {
						continue;
					}
					
					Next2_S = DesS[(X + PredDist) + (Y + 0) * DesPitch];
					Next2_V = DesV[(X + PredDist) + (Y + 0) * DesPitch];
					if (Y < SizeY1) {
						auto Temp = 0;
						if (Y == 0) {
							auto Curr0L_S = DesS[(X + PredDist) + (Y + 0) * DesPitch];
							auto Next1L_S = DesS[(X + PredDist) + (Y + 1) * DesPitch];
							auto Next2L_S = DesS[(X + PredDist) + (Y + 2) * DesPitch];
							Temp = DeltaQubicWaveletPre2(Curr0L_S, Next1L_S, Next2L_S);
						} else
						if (Y == 1) {
							auto Last1L_S = DesS[(X + PredDist) + (Y - 1) * DesPitch];
							auto Curr0L_S = DesS[(X + PredDist) + (Y + 0) * DesPitch];
							auto Next1L_S = DesS[(X + PredDist) + (Y + 1) * DesPitch];
							auto Next2L_S = DesS[(X + PredDist) + (Y + 2) * DesPitch];
							Temp = DeltaQubicWaveletPre1(Last1L_S, Curr0L_S, Next1L_S, Next2L_S);
						} else
						if (Y + PredDist < SizeY0) {
							auto Last2L_S = DesS[(X + PredDist) + (Y - 2) * DesPitch];
							auto Last1L_S = DesS[(X + PredDist) + (Y - 1) * DesPitch];
							auto Curr0L_S = DesS[(X + PredDist) + (Y + 0) * DesPitch];
							auto Next1L_S = DesS[(X + PredDist) + (Y + 1) * DesPitch];
							auto Next2L_S = DesS[(X + PredDist) + (Y + 2) * DesPitch];
							Temp = DeltaQubicWavelet(Last2L_S, Last1L_S, Curr0L_S, Next1L_S, Next2L_S);
						} else
						if (Y + 2 == SizeY0) {
							auto Last2L_S = DesS[(X + PredDist) + (Y - 2) * DesPitch];
							auto Last1L_S = DesS[(X + PredDist) + (Y - 1) * DesPitch];
							auto Curr0L_S = DesS[(X + PredDist) + (Y + 0) * DesPitch];
							auto Next1L_S = DesS[(X + PredDist) + (Y + 1) * DesPitch];
							Temp = DeltaQubicWaveletPost1(Last2L_S, Last1L_S, Curr0L_S, Next1L_S);
						} else
						if (Y + 1 == SizeY0) {
							auto Last2L_S = DesS[(X + PredDist) + (Y - 2) * DesPitch];
							auto Last1L_S = DesS[(X + PredDist) + (Y - 1) * DesPitch];
							auto Curr0L_S = DesS[(X + PredDist) + (Y + 0) * DesPitch];
							Temp = DeltaQubicWaveletPost2(Last2L_S, Last1L_S, Curr0L_S);
						}
						auto V = DesV[(X + PredDist) + Y * DesPitch] - Temp;
						DesV[(X + PredDist) + Y * DesPitch] = (tInt16)V;
					}
				}
				
				if (Y < 0 || X < 0 || (SizeX & SizeY) < 8) {
					continue;
				}
				
				auto TempH = 0;
				auto TempD = 0;
				if (X == 0) {
					TempH = DeltaQubicWaveletPre2(Curr0_S, Next1_S, Next2_S);
					TempD = DeltaQubicWaveletPre2(Curr0_V, Next1_V, Next2_V);
				} else
				if (X == 1) {
					TempH = DeltaQubicWaveletPre1(Last1_S, Curr0_S, Next1_S, Next2_S);
					TempD = DeltaQubicWaveletPre1(Last1_V, Curr0_V, Next1_V, Next2_V);
				} else
				if (X + PredDist < SizeX1) {
					TempH = DeltaQubicWavelet(Last2_S, Last1_S, Curr0_S, Next1_S, Next2_S);
					TempD = DeltaQubicWavelet(Last2_V, Last1_V, Curr0_V, Next1_V, Next2_V);
				} else
				if (X + 2 == SizeX1) {
					TempH = DeltaQubicWaveletPost1(Last2_S, Last1_S, Curr0_S, Next1_S);
					TempD = DeltaQubicWaveletPost1(Last2_V, Last1_V, Curr0_V, Next1_V);
				} else
				if (X + 1 == SizeX1) {
					TempH = DeltaQubicWaveletPost2(Last2_S, Last1_S, Curr0_S);
					TempD = DeltaQubicWaveletPost2(Last2_V, Last1_V, Curr0_V);
				}
				{
					auto H = DesH[X + Y * DesPitch] - TempH;
					DesH[X + Y * DesPitch] = (tInt16)H;
				}
				if (Y < SizeY1) {
					auto D = DesD[X + Y * DesPitch] - ((TempD + DFactor) >> DFactor);
					DesD[X + Y * DesPitch] = (tInt16)D;
				}
			}
		}
		
		STOP_CLOCK(SplitWithQubPred);
	}
	
	//=============================================================================
	static
	void ComposeWithQubPred(
		tLevel* SrcLevel,
		tLevel* DesLevel
	//=============================================================================
	) {
		START_CLOCK(ComposeWithQubPred);
		
		auto SizeX    = (tInt32)DesLevel->SizeXLeft;
		auto SizeY    = (tInt32)DesLevel->SizeYTop;
		auto DesPitch = (tInt32)DesLevel->Pitch;
		auto SrcPitch = (tInt32)SrcLevel->Pitch;
		auto Des      = DesLevel->S;
		auto SrcS     = SrcLevel->S;
		auto SrcH     = SrcLevel->H;
		auto SrcV     = SrcLevel->V;
		auto SrcD     = SrcLevel->D;
		
		auto SizeX0 = (tInt32)(SizeX + 1) >> 1;
		auto SizeX1 = (tInt32)(SizeX + 0) >> 1;
		auto SizeY0 = (tInt32)(SizeY + 1) >> 1;
		auto SizeY1 = (tInt32)(SizeY + 0) >> 1;
		
		auto Last2_S = 0;
		auto Last2_V = 0;
		auto Last1_S = 0;
		auto Last1_V = 0;
		auto Curr0_S = 0;
		auto Curr0_V = 0;
		auto Next1_S = 0;
		auto Next1_V = 0;
		auto Next2_S = 0;
		auto Next2_V = 0;
		
		for (auto Y = 0; Y < SizeY0; Y++) {
			auto Y2 = Y << 1;
			for (auto X = -PredDist; X < SizeX0; X++) {
				Last2_S = Last1_S;
				Last2_V = Last1_V;
				Last1_S = Curr0_S;
				Last1_V = Curr0_V;
				Curr0_S = Next1_S;
				Curr0_V = Next1_V;
				Next1_S = Next2_S;
				Next1_V = Next2_V;
				if (X + PredDist < SizeX0 && (SizeX & SizeY) >= 8) {
					Next2_S = SrcS[(tNat16)(X + PredDist) + Y * SrcPitch];
					
					if (Y < SizeY1) {
						auto Temp = 0;
						if (Y == 0) {
							auto Curr0L_S = SrcS[(X + PredDist) + (Y + 0) * SrcPitch];
							auto Next1L_S = SrcS[(X + PredDist) + (Y + 1) * SrcPitch];
							auto Next2L_S = SrcS[(X + PredDist) + (Y + 2) * SrcPitch];
							Temp = DeltaQubicWaveletPre2(Curr0L_S, Next1L_S, Next2L_S);
						} else
						if (Y == 1) {
							auto Last1L_S = SrcS[(X + PredDist) + (Y - 1) * SrcPitch];
							auto Curr0L_S = SrcS[(X + PredDist) + (Y + 0) * SrcPitch];
							auto Next1L_S = SrcS[(X + PredDist) + (Y + 1) * SrcPitch];
							auto Next2L_S = SrcS[(X + PredDist) + (Y + 2) * SrcPitch];
							Temp = DeltaQubicWaveletPre1(Last1L_S, Curr0L_S, Next1L_S, Next2L_S);
						} else
						if (Y + 2 < SizeY0) {
							auto Last2L_S = SrcS[(X + PredDist) + (Y - 2) * SrcPitch];
							auto Last1L_S = SrcS[(X + PredDist) + (Y - 1) * SrcPitch];
							auto Curr0L_S = SrcS[(X + PredDist) + (Y + 0) * SrcPitch];
							auto Next1L_S = SrcS[(X + PredDist) + (Y + 1) * SrcPitch];
							auto Next2L_S = SrcS[(X + PredDist) + (Y + 2) * SrcPitch];
							Temp = DeltaQubicWavelet(Last2L_S, Last1L_S, Curr0L_S, Next1L_S, Next2L_S);
						} else
						if (Y + 2 == SizeY0) {
							auto Last2L_S = SrcS[(X + PredDist) + (Y - 2) * SrcPitch];
							auto Last1L_S = SrcS[(X + PredDist) + (Y - 1) * SrcPitch];
							auto Curr0L_S = SrcS[(X + PredDist) + (Y + 0) * SrcPitch];
							auto Next1L_S = SrcS[(X + PredDist) + (Y + 1) * SrcPitch];
							Temp = DeltaQubicWaveletPost1(Last2L_S, Last1L_S, Curr0L_S, Next1L_S);
						} else
						if (Y + 1 == SizeY0) {
							auto Last2L_S = SrcS[(X + PredDist) + (Y - 2) * SrcPitch];
							auto Last1L_S = SrcS[(X + PredDist) + (Y - 1) * SrcPitch];
							auto Curr0L_S = SrcS[(X + PredDist) + (Y + 0) * SrcPitch];
							Temp = DeltaQubicWaveletPost2(Last2L_S, Last1L_S, Curr0L_S);
						}
						SrcV[(X + PredDist) + Y * SrcPitch] += (tInt16)Temp;
					}
					Next2_V = SrcV[(X + PredDist) + (Y + 0) * SrcPitch];
				}
				
				if (X >= 0 && (SizeX & SizeY) >= 8) {
					auto TempH = 0;
					auto TempD = 0;
					if (X == 0) {
						TempH = DeltaQubicWaveletPre2(Curr0_S, Next1_S, Next2_S);
						TempD = DeltaQubicWaveletPre2(Curr0_V, Next1_V, Next2_V);
					} else
					if (X == 1) {
						TempH = DeltaQubicWaveletPre1(Last1_S, Curr0_S, Next1_S, Next2_S);
						TempD = DeltaQubicWaveletPre1(Last1_V, Curr0_V, Next1_V, Next2_V);
					} else
					if (X + PredDist < SizeX1) {
						TempH = DeltaQubicWavelet(Last2_S, Last1_S, Curr0_S, Next1_S, Next2_S);
						TempD = DeltaQubicWavelet(Last2_V, Last1_V, Curr0_V, Next1_V, Next2_V);
					} else
					if (X + 2 == SizeX1) {
						TempH = DeltaQubicWaveletPost1(Last2_S, Last1_S, Curr0_S, Next1_S);
						TempD = DeltaQubicWaveletPost1(Last2_V, Last1_V, Curr0_V, Next1_V);
					} else
					if (X + 1 == SizeX1) {
						TempH = DeltaQubicWaveletPost2(Last2_S, Last1_S, Curr0_S);
						TempD = DeltaQubicWaveletPost2(Last2_V, Last1_V, Curr0_V);
					}
					SrcH[X + Y * SrcPitch] += (tInt16)TempH;
					if (Y < SizeY1) {
						SrcD[X + Y * SrcPitch] += (tInt16)(TempD + DFactor) >> DFactor;
					}
				}
				
				if (X >= 0) {
					auto X2 = X << 1;
					
					HaarWavelet::tQuad Arg;
					Arg.A = 0;
					Arg.B = 0;
					Arg.C = 0;
					Arg.D = 0;
					
					auto PosSrc = X + Y * SrcPitch;
					
					if (Y2 + 1 == SizeY) {
						if (X2 + 1 == SizeX) {
							Arg.A = SrcS[PosSrc];
							auto Res = HaarWavelet::DeCode(Arg);
							ASSERT(Res.A >= 0 &&
								Res.B >= 0 &&
								Res.C >= 0 &&
								Res.D >= 0
							);
							Des[(tNat32)((X2 + 0) + (Y2 + 0) * DesPitch)] = Res.A;
						} else {
							Arg.A = SrcS[PosSrc];
							Arg.B = SrcH[PosSrc];
							auto Res = HaarWavelet::DeCode(Arg);
							ASSERT(
								Res.A >= 0 &&
								Res.B >= 0 &&
								Res.C >= 0 &&
								Res.D >= 0
							);
							Des[(tNat32)((X2 + 0) + (Y2 + 0) * DesPitch)] = Res.A;
							Des[(tNat32)((X2 + 1) + (Y2 + 0) * DesPitch)] = Res.B;
						}
					} else {
						if (X2 + 1 == SizeX) {
							Arg.A = SrcS[PosSrc];
							Arg.C = SrcV[PosSrc];
							auto Res = HaarWavelet::DeCode(Arg);
							ASSERT(
								Res.A >= 0 &&
								Res.B >= 0 &&
								Res.C >= 0 &&
								Res.D >= 0
							);
							Des[(tNat32)((X2 + 0) + (Y2 + 0) * DesPitch)] = Res.A;
							Des[(tNat32)((X2 + 0) + (Y2 + 1) * DesPitch)] = Res.C;
						} else {
							Arg.A = SrcS[PosSrc];
							Arg.B = SrcH[PosSrc];
							Arg.C = SrcV[PosSrc];
							Arg.D = SrcD[PosSrc];
							auto Res = HaarWavelet::DeCode(Arg);
// TODO
//							ASSERT(
//								Res.A >= 0 &&
//								Res.B >= 0 &&
//								Res.C >= 0 &&
//								Res.D >= 0
//							);
							Des[(tNat32)((X2 + 0) + (Y2 + 0) * DesPitch)] = Res.A;
							Des[(tNat32)((X2 + 1) + (Y2 + 0) * DesPitch)] = Res.B;
							Des[(tNat32)((X2 + 0) + (Y2 + 1) * DesPitch)] = Res.C;
							Des[(tNat32)((X2 + 1) + (Y2 + 1) * DesPitch)] = Res.D;
						}
					}
				}
			}
		}
		
		STOP_CLOCK(ComposeWithQubPred);
	}
	
#	ifdef TEST_
		//=============================================================================
		static inline
		void Test(
		//=============================================================================
		) {
			const tInt32 SizeMax = 128;
			const tInt32 DeltaPitch = 16;
			
			tInt16 Big[(SizeMax + DeltaPitch) * SizeMax] = {};
			tInt16 S[(SizeMax/2 + DeltaPitch) * SizeMax/2] = {};
			tInt16 H[(SizeMax/2 + DeltaPitch) * SizeMax/2] = {};
			tInt16 V[(SizeMax/2 + DeltaPitch) * SizeMax/2] = {};
			tInt16 D[(SizeMax/2 + DeltaPitch) * SizeMax/2] = {};
			
			for (auto I = 1; I < 500; I += 1) {
				for (auto Y = 0; Y < SizeMax; Y += 1) {
					for (auto X = 0; X < SizeMax; X += 1) {
						Big[X + Y*SizeMax] = ((X+257) * (Y+263) * I) % 251;
					}
				}
				
				for (auto Size = 2; Size <= SizeMax; Size <<= 1) {
					SplitWithQubPred  (Size, Size, Size+DeltaPitch, Size/2+DeltaPitch, Big, S, H, V, D);
					ComposeWithQubPred(Size, Size, Size+DeltaPitch, Size/2+DeltaPitch, Big, S, H, V, D);
					
					for (auto Y = 0; Y < SizeMax; Y += 1) {
						for (auto X = 0; X < SizeMax; X += 1) {
							auto A = Big[X + Y*SizeMax];
							auto B = ((X+257) * (Y+263) * I) % 251;
							ASSERT(A == B);
						}
					}
				}
			}
		}
#	endif
}

namespace BmpHelper {
	
	using tFileHeader = struct {
		tNat32 Size;   // size of file
		tNat32 _;      // reserved
		tNat32 Offset; // sizeof(MagicString) + sizeof(tBmpFileHeader) + sizeof(tBmpInfoHeader) + sizeof(mask) + sizeof(palette)
	};
	
	using tInfoHeader = struct {
		tNat32 HeaderSize;     // 40
		tInt32 SizeX;
		tInt32 SizeY;          // negativ -> top-down; positiv -> bottom-up
		tNat16 _;              // reserved
		tNat16 BitCount;       // Pallette: 1, 4, 8; RGB: 16, 24, 32
		tNat32 Compression;    // 0 -> off; 3 -> Masken (bei 16 & 32 bpp); 1|2 -> Lauflangenkodierung bei 8bpp|4bpp
		tNat32 SizeImage;      // if Compression == 1 or 2
		tNat32 PixelPerMeterX;
		tNat32 PixelPerMeterY;
		tNat32 ClrUsed;        // Farben in Farbtabelle (0=max)
		tNat32 ClrImportant;   // Anzahl verwendeter Farben aus der Farbtabelle (0=max)
	};
	
	//=============================================================================
	static inline
	void GetColorMask(
		tStream* StreamIn,
		tNat32*  MaskR,
		tNat32*  MaskG,
		tNat32*  MaskB
	//=============================================================================
	) {
		if (
			!fread(&MaskR, sizeof(tNat32), 1, StreamIn) ||
			!fread(&MaskG, sizeof(tNat32), 1, StreamIn) ||
			!fread(&MaskB, sizeof(tNat32), 1, StreamIn)
		) {
			fprintf(gLogStream, "ERROR: can't read file");
			exit(-1);
		}
	}
	
	//=============================================================================
	static inline
	tNat32 GetShift(
		tNat32 Mask
	//=============================================================================
	) {
		for (auto Shift = 0u; Shift < 32; Shift += 1) {
			if (Mask & 1) {
				return Shift;
			}
			Mask >>= 1;
		}
		return 0;
	};
	
	//=============================================================================
	static
	void GetLayerFromBmp16(
		tStream*       StreamIn,
		tInfoHeader    BmpInfoHeader,
		Layer::tLevel* LayerA,
		Layer::tLevel* LayerY,
		Layer::tLevel* LayerCg,
		Layer::tLevel* LayerCo
	//=============================================================================
	) {
		START_CLOCK(GetLayerFromBmp16)
		
		auto MaskR = 0x00007C00u;
		auto MaskG = 0x000003E0u;
		auto MaskB = 0x0000001Fu;
		if (BmpInfoHeader.Compression == 3) {
			GetColorMask(StreamIn, &MaskR, &MaskG, &MaskB);
		}
		auto MaskA = ~(MaskR | MaskG | MaskB) & 0x0000FFFFu;
		
		auto ShiftA = GetShift(MaskA);
		auto ShiftR = GetShift(MaskR);
		auto ShiftG = GetShift(MaskG);
		auto ShiftB = GetShift(MaskB);
		
		auto SizeX = BmpInfoHeader.SizeX;
		auto RowSize = (tSize)(2*SizeX + 3) & (~3);
		auto Row = (tNat16*)malloc(RowSize); // TODO: use FixSize row to skip malloc
		for (auto Y = 0; fread(Row, RowSize, 1, StreamIn); Y += 1) {
			for (auto X = 0; X < SizeX; X += 1) {
				auto J = Y*SizeX + X;
				LayerA->S[J] = (tInt16)(Row[X] & MaskA) >> ShiftA;
				
				auto R = (tInt16)(Row[X] & MaskR) >> ShiftR;
				auto G = (tInt16)(Row[X] & MaskG) >> ShiftG;
				auto B = (tInt16)(Row[X] & MaskB) >> ShiftB;
				
				auto RB = R + B;
				
				LayerY->S[J]  = (tInt16)((G<<1) + RB + 2) >> 2;
				LayerCg->S[J] = (tInt16)((G<<1) - RB + 1) >> 1;
				LayerCo->S[J] = (tInt16)(B - R);
			}
		}
		free(Row);
		STOP_CLOCK(GetLayerFromBmp16)
	}
	
	//=============================================================================
	static
	void GetLayerFromBmp24(
		tStream*       StreamIn,
		tInfoHeader    BmpInfoHeader,
		Layer::tLevel* LayerA,
		Layer::tLevel* LayerY,
		Layer::tLevel* LayerCg,
		Layer::tLevel* LayerCo
	//=============================================================================
	) {
		START_CLOCK(GetLayerFromBmp24);
		auto SizeX = BmpInfoHeader.SizeX;
		auto Pitch = LayerA->Pitch;
		auto RowSize = (tSize)((3*SizeX + 3) & (~3));
		auto Row = (tNat8*)malloc(RowSize); // TODO: use FixSize row to skip malloc
		for (auto Y = 0; fread(Row, RowSize, 1, StreamIn); Y += 1) {
			auto I = 0;
			for (auto X = 0; X < SizeX; X += 1) {
				auto J = Y*Pitch + X;
				LayerA->S[J] = 255;
				
				auto B = Row[I++];
				auto G = Row[I++];
				auto R = Row[I++];
				
				auto RB = R + B;
				
				LayerY->S[J]  = (tInt16)((G<<1) + RB + 2) >> 2;
				LayerCg->S[J] = (tInt16)((G<<1) - RB + 1) >> 1;
				LayerCo->S[J] = (tInt16)(B - R);
			}
		}
		free(Row);
		
		STOP_CLOCK(GetLayerFromBmp24);
	}
	
	//=============================================================================
	static
	void GetLayerFromBmp32(
		tStream*       StreamIn,
		tInfoHeader    BmpInfoHeader,
		Layer::tLevel* LayerA,
		Layer::tLevel* LayerY,
		Layer::tLevel* LayerCg,
		Layer::tLevel* LayerCo
	//=============================================================================
	) {
		START_CLOCK(GetLayerFromBmp32);
		
		auto MaskR = 0x00FF0000u;
		auto MaskG = 0x0000FF00u;
		auto MaskB = 0x000000FFu;
		if (BmpInfoHeader.Compression == 3) {
			GetColorMask(StreamIn, &MaskR, &MaskG, &MaskB);
		}
		auto MaskA = ~(MaskR | MaskG | MaskB);
		
		auto ShiftA = GetShift(MaskA);
		auto ShiftR = GetShift(MaskR);
		auto ShiftG = GetShift(MaskG);
		auto ShiftB = GetShift(MaskB);
		
		auto SizeX = BmpInfoHeader.SizeX;
		auto RowSize = (tSize)(4*SizeX);
		auto Row = (tNat32*)malloc(RowSize); // TODO: use FixSize row to skip malloc
		for (auto Y = 0; fread(Row, RowSize, 1, StreamIn); Y += 1) {
			for (auto X = 0; X < SizeX; X += 1) {
				auto J = Y*SizeX + X;
				LayerA->S[J] = (tInt16)(Row[X] & MaskA) >> ShiftA;
				
				auto R = (tInt16)(Row[X] & MaskR) >> ShiftR;
				auto G = (tInt16)(Row[X] & MaskG) >> ShiftG;
				auto B = (tInt16)(Row[X] & MaskB) >> ShiftB;
				
				auto RB = R + B;
				
				LayerY->S[J]  = (tInt16)((G<<1) + RB + 2) >> 2;
				LayerCg->S[J] = (tInt16)((G<<1) - RB + 1) >> 1;
				LayerCo->S[J] = (tInt16)(B - R);
			}
		}
		free(Row);
		
		STOP_CLOCK(GetLayerFromBmp32);
	}
	
	//=============================================================================
	static
	void PutLayerToBmp32(
		tStream*       StreamOut,
		tInfoHeader    BmpInfoHeader,
		Layer::tLevel* LayerA,
		Layer::tLevel* LayerY,
		Layer::tLevel* LayerCg,
		Layer::tLevel* LayerCo
	//=============================================================================
	) {
		START_CLOCK(PutLayerToBmp32);
		
		tNat32 Mask[3];
		auto MaskR = Mask[0] = 0x00FF0000u; // R
		auto MaskG = Mask[1] = 0x0000FF00u; // G
		auto MaskB = Mask[2] = 0x000000FFu; // B
		auto MaskA = ~(MaskR | MaskG | MaskB);
		if (!fwrite(&Mask, sizeof(Mask), 1, StreamOut)) {
			fprintf(gLogStream, "ERROR: fail writing to output!!!");
			exit(-1);
		}
		
		auto ShiftA = GetShift(MaskA);
		auto ShiftR = GetShift(MaskR);
		auto ShiftG = GetShift(MaskG);
		auto ShiftB = GetShift(MaskB);
		
		auto SizeX = BmpInfoHeader.SizeX;
		auto Pitch = LayerA->Pitch;
		auto RowSize = (tSize)(4*SizeX);
		auto Row = (tNat32*)malloc(RowSize); // TODO: use FixSize row to skip malloc
		for (auto Y = 0; Y < BmpInfoHeader.SizeY; Y += 1) {
			for (auto X = 0; X < SizeX; X += 1) {
				auto J = Y*Pitch + X;
				
				auto A  = LayerA->S[J];
				auto Y_  = LayerY->S[J];
				auto Cg = LayerCg->S[J];
				auto Co = LayerCo->S[J];
				
				auto RB = (Y_<<1) - Cg;
				
				auto G = ((Y_<<1) + Cg) >> 1;
				auto B = (RB + Co + 1) >> 1;
				auto R = (RB - Co + 1) >> 1;
				
				auto Value = (
					Min(Max(0, R), 255) << ShiftR |
					Min(Max(0, G), 255) << ShiftG |
					Min(Max(0, B), 255) << ShiftB |
					Min(Max(0, A), 255) << ShiftA
				);
				Row[X] = Value;
			}
			fwrite(Row, RowSize, 1, StreamOut);
		}
		free(Row);
		
		STOP_CLOCK(PutLayerToBmp32);
	}
}

using iCurve = struct {
	void* Env;
	tFunc3<void*, tNat32, tNat32, void>* Init;
	tFunc1<void*, void>* Next;
	tFunc1<void*, tNat32>* GetX;
	tFunc1<void*, tNat32>* GetY;
};

//=============================================================================
static inline
void Init(
	iCurve* CurveInterface,
	tNat32  SizeX,
	tNat32  SizeY
//=============================================================================
) {
	CurveInterface->Init(CurveInterface->Env, SizeX, SizeY);
}

//=============================================================================
static inline
void Next(
	iCurve* CurveInterface
//=============================================================================
) {
	CurveInterface->Next(CurveInterface->Env);
}

//=============================================================================
static inline
tNat32 GetX(
	iCurve* CurveInterface
//=============================================================================
) {
	return CurveInterface->GetX(CurveInterface->Env);
}

//=============================================================================
static inline
tNat32 GetY(
	iCurve* CurveInterface
//=============================================================================
) {
	return CurveInterface->GetY(CurveInterface->Env);
}

namespace Scanline {
	using tState = struct {
		tNat32 MaxX;
		tNat32 X;
		tNat32 Y;
		tNat32 _;
	};
	
	//=============================================================================
	static inline
	void Init(
		tState* State,
		tNat32  SizeX,
		tNat32  SizeY
	//=============================================================================
	) {
		*State = {};
		State->MaxX = SizeX;
	}
	
	//=============================================================================
	static inline
	void Next(
		tState* State
	//=============================================================================
	) {
		State->X += 1;
		if (State->X > State->MaxX) {
			State->Y += 1;
			State->X  = 0;
		}
	}
	
	//=============================================================================
	static inline
	tNat32 GetX(
		tState* State
	//=============================================================================
	) {
		return State->X;
	}
	
	//=============================================================================
	static inline
	tNat32 GetY(
		tState* State
	//=============================================================================
	) {
		return State->Y;
	}
	
	//=============================================================================
	static inline
	iCurve GetInterface(
		tState* State
	//=============================================================================
	) {
		iCurve Interface = {};
		
		Interface.Env = State;
		Interface.Init = (tFunc3<void*, tNat32, tNat32, void>*)Init;
		Interface.Next = (tFunc1<void*, void>*)Next;
		Interface.GetX = (tFunc1<void*, tNat32>*)GetX;
		Interface.GetY = (tFunc1<void*, tNat32>*)GetY;
		
		return Interface;
	}
}

namespace ZCurve {
	using tState = struct {
		tNat64 Step;
		tNat32 X;
		tNat32 Y;
	};
	
	//=============================================================================
	static inline
	void Init(
		tState* State,
		tNat32  SizeX,
		tNat32  SizeY
	//=============================================================================
	) {
		*State = {};
	}
	
	//=============================================================================
	static inline
	void Next(
		tState* State
	//=============================================================================
	) {
		State->Step += 1;
		auto Step = State->Step;
		
		auto XY = (Step & 0x55555555) | ((Step & 0xAAAAAAAA) << 31);
		XY = (XY & 0x1111111111111111) | ((XY & 0x4444444444444444) >> 1);
		XY = (XY & 0x0303030303030303) | ((XY & 0x3030303030303030) >> 2);
		XY = (XY & 0x000F000F000F000F) | ((XY & 0x0F000F000F000F00) >> 4);
		XY = (XY & 0x000000FF000000FF) | ((XY & 0x00FF000000FF0000) >> 4);
		
		State->X = XY & 0xFFFF;
		State->Y = (XY >> 32) & 0xFFFF;
	}
	
	//=============================================================================
	static inline
	tNat32 GetX(
		tState* State
	//=============================================================================
	) {
		return State->X;
	}
	
	//=============================================================================
	static inline
	tNat32 GetY(
		tState* State
	//=============================================================================
	) {
		return State->Y;
	}
	
	//=============================================================================
	static inline
	iCurve GetInterface(
		tState* State
	//=============================================================================
	) {
		iCurve Interface = {};
		
		Interface.Env = State;
		Interface.Init = (tFunc3<void*, tNat32, tNat32, void>*)Init;
		Interface.Next = (tFunc1<void*, void>*)Next;
		Interface.GetX = (tFunc1<void*, tNat32>*)GetX;
		Interface.GetY = (tFunc1<void*, tNat32>*)GetY;
		
		return Interface;
	}
	
}

namespace HilbertCurve {
	
	using tStack = tNat64;
	
	using tState = struct {
		tStack Stack;
		tNat32 X;
		tNat32 Y;
	};
	
	using tStepType = enum {
		A = 0,
		B = 1,
		C = 2,
		D = 3,
	};
	
	using tStepNumber = enum {
		Alpha = 3,
		Betha = 2,
		Gamma = 1,
		End = 0,
	};
	
	using tDirection = enum {
		N = 0,
		W = 1,
		S = 2,
		E = 3,
	};
	
	//=============================================================================
	inline static
	tStepType StepType(
		tStack Stack
	//=============================================================================
	) {
		return (tStepType)(Stack & 3);
	}
	
	//=============================================================================
	inline static
	tStepNumber StepNumber(
		tStack Stack
	//=============================================================================
	) {
		return (tStepNumber)((Stack >> 2) & 3);
	}
	
	//=============================================================================
	inline static
	void AddLevel(
		tStack* Stack
	//=============================================================================
	) {
		static tNat32 Matrix[] = { D, A, A, B };
		
		auto OldStepType = StepType(*Stack);
		auto OldStepNumber = StepNumber(*Stack);
		
		auto NewStepType = Matrix[OldStepNumber] ^ OldStepType;
		const tNat32 NewStepNumber = Alpha;
		
		*Stack <<= 4;
		*Stack |= (NewStepNumber << 2) | NewStepType;
	}
	
	//=============================================================================
	inline static
	tStack New(
		tNat32 Levels
	//=============================================================================
	) {
		auto Stack = (tStack)((Alpha << 2) | A);
		
		while (Levels --> 0) {
			Stack = (Stack << 4) | (Alpha << 2) | ((Stack & 0x3) ^ 1);
		}
		
		return Stack;
	}
	
	//=============================================================================
	inline static
	void Init(
		tState* State,
		tNat32  SizeX,
		tNat32  SizeY
	//=============================================================================
	) {
		auto Levels = 0;
		for (auto Temp = SizeX | SizeY; Temp; Temp >>= 1) {
			Levels += 1;
		}
		
		auto Stack = (tStack)((Alpha << 2) | A);
		while (Levels --> 0) {
			Stack = (Stack << 4) | (Alpha << 2) | ((Stack & 0x3) ^ 1);
		}
		
		State->Stack = Stack;
		State->X = 0;
		State->Y = 0;
	}
	
	//=============================================================================
	inline static
	tDirection Direction(
		tStack Stack
	//=============================================================================
	) {
		return (tDirection)((Stack ^ (Stack >> 2)) & 3);
	}
	
	//=============================================================================
	inline static
	void Next(
		tState* State
	//=============================================================================
	) {
//		START_CLOCK(HilbertCurve_Next);
		
		auto RemovedLevels = 0;
		while (StepNumber(State->Stack) == End) {
			State->Stack >>= 4; // RemoveLevel
			RemovedLevels += 1;
		}
		
		auto Result = Direction(State->Stack);
		State->Stack -= 1 << 2; // NextStepNumber
		
		while (RemovedLevels --> 0) { AddLevel(&State->Stack); }
		
		switch (Result) {
			case N: {
				State->Y -= 1;
			} break;
			case S: {
				State->Y += 1;
			} break;
			case E: {
				State->X += 1;
			} break;
			case W: {
				State->X -= 1;
			} break;
			default: {
				ASSERT(false);
			}
		}
		
//		STOP_CLOCK(HilbertCurve_Next);
	}
	
	//=============================================================================
	inline static
	tNat32 GetX(
		tState* State
	//=============================================================================
	) {
		return State->X;
	}
	
	//=============================================================================
	inline static
	tNat32 GetY(
		tState* State
	//=============================================================================
	) {
		return State->Y;
	}
	
	//=============================================================================
	inline static
	iCurve GetInterface(
		tState* State
	//=============================================================================
	) {
		iCurve Interface;
		
		Interface.Env  = State;
		Interface.Init = (tFunc3<void*, tNat32, tNat32, void>*)Init;
		Interface.Next = (tFunc1<void*, void>*)Next;
		Interface.GetX = (tFunc1<void*, tNat32>*)GetX;
		Interface.GetY = (tFunc1<void*, tNat32>*)GetY;
		
		return Interface;
	}
	
#	ifdef TEST
		//=============================================================================
		static inline
		void Test(
		//=============================================================================
		) {
			const tNat32 Bits = 10;
			
			auto Array = (tBool*)malloc(1 << (2*Bits));
			for (auto I = 1 << (2*Bits); I --> 0;) {
				Array[I] = false;
			}
			Array[0] = true;
			
			tState HilbertState = {};
			HilbertState.Stack = HilbertCurve::New(Bits);
			
			Array[0] = true;
			
			for (auto Steps = (1 << (2*Bits)) - 1; Steps --> 0;) {
				HilbertCurve::Next(&HilbertState);
				ASSERT(Array[HilbertState.X | (HilbertState.Y << Bits)] == false);
				Array[HilbertState.X | (HilbertState.Y << Bits)] = true;
			}
			
			free(Array);
		}
#	endif
}

#if 1
#	define ReadNat(BitReader)       FibonacciCode::Read(BitReader)
#	define WriteNat(BitWriter, Nat) FibonacciCode::Write((BitWriter), (Nat))
#elif 0
#	define ReadNat(BitReader)       EliasGammaCode::Read(BitReader)
#	define WriteNat(BitWriter, Nat) EliasGammaCode::Write((BitWriter), (Nat))
#else
#	define ReadNat(BitReader)       ReadNat_(BitReader)
#	define WriteNat(BitWriter, Nat) WriteNat_((BitWriter), (Nat))
#endif

#if 1
	using tCurveState = HilbertCurve::tState;
#elif 1
	using tCurveState = ZCurve::tState;
#else
	using tCurveState = Scanline::tState;
#endif

//=============================================================================
static
void WriteLayer(
	BufferdStream::tWriter* Stream,
	tInt16*                 DataPtr,
	tNat32                  SizeX,
	tNat32                  SizeY,
	tNat32                  Pitch,
	tNat32                  Compression
//=============================================================================
) {
	START_CLOCK(WriteLayer);
	
	BitStream::tWriter BitWriter;
	BitStream::Init(&BitWriter, Stream);
	
	auto BitWriterInterface = GetInterface(&BitWriter);
	
//	ArithmeticBitStream::tWriter ABitWriter;
//	ArithmeticBitStream::Init(&ABitWriter, &BitWriterInterface_);
//	
//	auto BitWriterInterface = GetInterface(&ABitWriter);
	
	WriteNat(&BitWriterInterface, Compression);
	
	auto Levels = (tNat32)0;
	for (auto Size = Max(SizeX, SizeY); Size > 0; Size >>= 1) {
		Levels += 1;
	}
	
	auto Zeros = (tNat32)0;
	tCurveState CurveState;
	auto Curve = GetInterface(&CurveState);
	Init(&Curve, SizeX, SizeY);
	
	auto CompressionOffset = (1 << Compression) >> 1;
	
	tNat32 HistValues[1<<15] = {};
	tNat32 HistZeros[16] = {};
	
	for (auto I = (tNat32)0; I < SizeX*SizeY;) {
		auto X = GetX(&Curve);
		auto Y = GetY(&Curve);
		if (X < SizeX && Y < SizeY) {
			I += 1;
			auto Value = (DataPtr[X + Y*Pitch] + CompressionOffset) >> Compression;
			if (Value == 0) {
				Zeros += 1;
			} else {
				if (Zeros) {
					HistValues[0] += Zeros;
					HistZeros[UsedBits(Zeros)] += 1;
				}
				HistValues[Abs(Value)] += 1;
				
				if (Zeros) {
					if (Zeros == 1) {
						WriteNat(&BitWriterInterface, 0);
					} else {
						WriteNat(&BitWriterInterface, 1);
						WriteNat(&BitWriterInterface, Zeros - 2);
					}
					Zeros = 0;
					WriteNat(&BitWriterInterface, ((Abs(Value) - 1) << 1) + ((Value < 0) ? 1 : 0));
				} else {
					WriteNat(&BitWriterInterface, (Abs(Value) << 1) + ((Value < 0) ? 1 : 0));
				}
			}
		}
		
		Next(&Curve);
	}
	
	if (Zeros) {
		HistValues[0] += Zeros;
		HistZeros[UsedBits(Zeros)] += 1;
		if (Zeros == 1) {
			WriteNat(&BitWriterInterface, 0);
		} else {
			WriteNat(&BitWriterInterface, 1);
			WriteNat(&BitWriterInterface, Zeros - 2);
		}
	}
	
	Close(&BitWriterInterface);
	
	fprintf(gLogStream, "\n\n");
	fprintf(gLogStream, "%d x %d\n", SizeX, SizeY);
	for (auto Value = 0; Value < 1<<15; Value += 1) {
		auto Count = HistValues[Value];
		if (Count > 0) {
			fprintf(gLogStream, "V %d %d\n", Value, Count);
		}
	}
	fprintf(gLogStream, "--\n");
	for (auto Bits = 0; Bits < 16; Bits += 1) {
		auto Count = HistZeros[Bits];
		if (Count > 0) {
			fprintf(gLogStream, "Z %d %d\n", Bits, Count);
		}
	}
	
	STOP_CLOCK(WriteLayer);
}

//=============================================================================
static
void ReadLayer(
	BufferdStream::tReader* Stream,
	tInt16*                 DataPtr,
	tNat32                  SizeX,
	tNat32                  SizeY,
	tNat32                  Pitch
//=============================================================================
) {
	START_CLOCK(ReadLayer);
	
	BitStream::tReader BitStream;
	BitStream::Init(&BitStream, Stream);
	
	auto BitReaderInterface = GetInterface(&BitStream);
	
//	ArithmeticBitStream::tReader ABitReader;
//	ArithmeticBitStream::Init(&ABitReader, &BitReaderInterface_);
//	
//	auto BitReaderInterface = GetInterface(&ABitReader);
	
	auto Compression = (tNat32)ReadNat(&BitReaderInterface);
	
	auto Levels = (tNat32)0;
	for (auto Size = Max(SizeX, SizeY); Size > 0; Size >>= 1) {
		Levels += 1;
	}
	auto Zeros = (tNat64)0;
	auto WasLastZero = false;
	tCurveState CurveState;
	auto Curve = GetInterface(&CurveState);
	Init(&Curve, SizeX, SizeY);
	
	for (auto I = (tNat32)0; I < SizeX * SizeY;) {
		auto X = GetX(&Curve);
		auto Y = GetY(&Curve);
		if (X < SizeX && Y < SizeY) {
			I += 1;
			tNat32 Sign;
			tNat64 Value;
			if (Zeros) {
				Zeros -= 1;
				Sign = (tNat32)1;
				Value = (tNat64)0;
			} else if (WasLastZero) {
				Value = ReadNat(&BitReaderInterface);
				Sign = (Value & 1) ? -1 : 1;
				Value >>= 1;
				Value += 1;
			} else {
				Value = ReadNat(&BitReaderInterface);
				Sign = (Value & 1) ? -1 : 1;
				Value >>= 1;
				if (Value == 0 && Sign == -1) {
					Zeros = ReadNat(&BitReaderInterface) + 1;
				}
			}
			WasLastZero = (Value == 0);
			DataPtr[X + Y*Pitch] = (tInt16)(Sign * (Value << Compression));
		}
		Next(&Curve);
	}
	
	Close(&BitReaderInterface);
	
	STOP_CLOCK(ReadLayer);
}

//=============================================================================
void EnCode(
	tStream*                StreamIn,
	BufferdStream::tWriter* StreamOut,
	tNat32                  Compression
//=============================================================================
) {
	START_CLOCK(EnCode);
	
	char MagicString[2];
	if (!fread(&MagicString, sizeof(MagicString), 1, StreamIn)) {
		fprintf(gLogStream, "ERROR: can't read File");
		exit(-1);
	}
	if (
		MagicString[0] != 'B' ||
		MagicString[1] != 'M'
	) {
		fprintf(gLogStream, "ERROR: wrong MagicString %.2s", MagicString);
		exit(-1);
	}
	
	BmpHelper::tFileHeader BmpFileHeader;
	if (!fread(&BmpFileHeader, sizeof(BmpHelper::tFileHeader), 1, StreamIn)) {
		fprintf(gLogStream, "ERROR: can't read file");
		exit(-1);
	}
	
	BmpHelper::tInfoHeader BmpInfoHeader;
	if (!fread(&BmpInfoHeader, sizeof(BmpHelper::tInfoHeader), 1, StreamIn)) {
		fprintf(gLogStream, "ERROR: can't read file");
		exit(-1);
	}
	
	if (BmpInfoHeader.HeaderSize != 40) {
		fprintf(gLogStream, "Invalid fileformat");
		exit(-1);
	}
	
	auto SizeX = (tNat16)BmpInfoHeader.SizeX;
	auto SizeY = (tNat16)BmpInfoHeader.SizeY;
	
	Layer::tLevel Layers[4][32];
	auto LayerCount = 4u;
	auto BufferSize = 4 * SizeX * SizeY;
	auto Buffer = (tInt16*)malloc(LayerCount * BufferSize * sizeof(tInt16));
	
	auto InitLevel = 0u;
	{
		auto TempBuffer = Buffer;
		for (auto Layer = LayerCount; Layer --> 0;) {
			InitLevel = Layer::Init(&Layers[Layer][0], SizeX, SizeY, TempBuffer, BufferSize);
			TempBuffer += BufferSize;
		}
	}
	
	if (BmpInfoHeader.BitCount == 16) {
		BmpHelper::GetLayerFromBmp16(StreamIn, BmpInfoHeader, &Layers[0][InitLevel], &Layers[1][InitLevel], &Layers[2][InitLevel], &Layers[3][InitLevel]);
	} else
	if (BmpInfoHeader.BitCount == 24) {
		BmpHelper::GetLayerFromBmp24(StreamIn, BmpInfoHeader, &Layers[0][InitLevel], &Layers[1][InitLevel], &Layers[2][InitLevel], &Layers[3][InitLevel]);
	} else
	if (BmpInfoHeader.BitCount == 32) {
		BmpHelper::GetLayerFromBmp32(StreamIn, BmpInfoHeader, &Layers[0][InitLevel], &Layers[1][InitLevel], &Layers[2][InitLevel], &Layers[3][InitLevel]);
	} else {
		fprintf(gLogStream, "Invalid fileformat");
		exit(-1);
	}
	
	// Berechnen
	for (auto Layer = LayerCount; Layer --> 0;) {
		for (auto Level = InitLevel; Level --> 0;) {
			Layer::SplitWithQubPred(&Layers[Layer][Level+1], &Layers[Layer][Level]);
		}
	}
	
	{ // Output
		BitStream::tWriter BitWriter;
		BitStream::Init(&BitWriter, StreamOut);
		auto BitWriterInterface = GetInterface(&BitWriter);
		
		WriteNat(&BitWriterInterface, SizeX);
		WriteNat(&BitWriterInterface, SizeY);
		
		for (auto Layer = LayerCount; Layer --> 0;) {
			auto Value = *Layers[Layer][0].S;
			WriteBit(&BitWriterInterface, Value < 0);
			WriteNat(&BitWriterInterface, Abs(Value));
			
			ASSERT(Layers[Layer][0].SizeXLeft == 1);
			ASSERT(Layers[Layer][0].SizeYTop == 1);
		}
		
		Close(&BitWriterInterface);
		
		for (auto Level = (tNat32)0; Level < InitLevel; Level += 1) {
			auto SizeXLeft = Layers[0][Level].SizeXLeft;
			auto SizeXRight = Layers[0][Level].SizeXRight;
			auto SizeYTop = Layers[0][Level].SizeYTop;
			auto SizeYBottom = Layers[0][Level].SizeYBottom;
			auto Pitch = Layers[0][Level].Pitch;
			
			ASSERT(SizeXLeft >= SizeXRight);
			ASSERT(SizeYTop >= SizeYBottom);
			
			ASSERT(SizeXLeft + SizeXRight == Layers[0][Level+1].SizeXLeft);
			ASSERT(SizeYTop + SizeYBottom == Layers[0][Level+1].SizeYTop);
			
			for (auto Layer = (tNat32)0; Layer < LayerCount; Layer += 1) {
				auto LayerTemp = Layers[Layer][Level];
				
				auto Compression_ = Max(0, Level - InitLevel + Compression - 2*(Layer == 1));
				
				WriteLayer(StreamOut, LayerTemp.H, SizeXRight, SizeYTop,    Pitch, Compression_);
				WriteLayer(StreamOut, LayerTemp.V, SizeXLeft,  SizeYBottom, Pitch, Compression_);
				WriteLayer(StreamOut, LayerTemp.D, SizeXRight, SizeYBottom, Pitch, Compression_);
			}
		}
	}
	
	STOP_CLOCK(EnCode);
}

//=============================================================================
void DeCode(
	BufferdStream::tReader* StreamIn,
	tStream*                StreamOut
//=============================================================================
) {
	START_CLOCK(DeCode);
	
	BitStream::tReader BitReader;
	BitStream::Init(&BitReader, StreamIn);
	auto BitReaderInterface = GetInterface(&BitReader);
	
	auto SizeX = (tNat32)ReadNat(&BitReaderInterface);
	auto SizeY = (tNat32)ReadNat(&BitReaderInterface);
	
	Layer::tLevel Layers[4][32];
	auto LayerCount = 4u;
	auto InitLevel = 0;
	{
		auto BufferSize = 4 * SizeX * Max(SizeY, 8);
		auto Buffer = (tInt16*)malloc(LayerCount * BufferSize * sizeof(tInt16));
		
		{
			auto TempBuffer = Buffer;
			for (auto Layer = LayerCount; Layer --> 0;) {
				InitLevel = Layer::Init(&Layers[Layer][0], SizeX, SizeY, TempBuffer, BufferSize);
				TempBuffer += BufferSize;
			}
		}
		
		for (auto Layer = LayerCount; Layer --> 0;) {
			auto Sign = ReadBit(&BitReaderInterface) ? -1 : 1;
			*Layers[Layer][0].S = (tInt16)(Sign * ReadNat(&BitReaderInterface));
			
			ASSERT(Layers[Layer][0].SizeXLeft == 1);
			ASSERT(Layers[Layer][0].SizeYTop == 1);
		}
		
		Close(&BitReaderInterface);
		
		for (auto Level = 0; Level < InitLevel; Level += 1) {
			auto SizeXLeft = Layers[0][Level].SizeXLeft;
			auto SizeXRight = Layers[0][Level].SizeXRight;
			auto SizeYTop = Layers[0][Level].SizeYTop;
			auto SizeYBottom = Layers[0][Level].SizeYBottom;
			auto Pitch = Layers[0][Level].Pitch;
			
			ASSERT(SizeXLeft >= SizeXRight);
			ASSERT(SizeYTop >= SizeYBottom);
			
			ASSERT(SizeXLeft + SizeXRight == Layers[0][Level+1].SizeXLeft);
			ASSERT(SizeYTop + SizeYBottom == Layers[0][Level+1].SizeYTop);
			
			for (auto Layer = (tNat32)0; Layer < LayerCount; Layer += 1) {
				auto LayerTemp = Layers[Layer][Level];
				
				ReadLayer(StreamIn, LayerTemp.H, SizeXRight, SizeYTop,    Pitch);
				ReadLayer(StreamIn, LayerTemp.V, SizeXLeft,  SizeYBottom, Pitch);
				ReadLayer(StreamIn, LayerTemp.D, SizeXRight, SizeYBottom, Pitch);
			}
		}
	}
	
	{ // Berechnen
		for (auto Layer = LayerCount; Layer --> 0;) {
			for (auto Level = 0; Level < InitLevel; Level += 1) {
				Layer::ComposeWithQubPred(&Layers[Layer][Level], &Layers[Layer][Level+1]);
			}
		}
	}
	
	{ // Output
		tChar8 MagicString[2];
		MagicString[0] = 'B';
		MagicString[1] = 'M';
		if (!fwrite(&MagicString, sizeof(MagicString), 1, StreamOut)) {
			fprintf(gLogStream, "ERROR: fail writing to output!!!");
			exit(-1);
		}
		
		BmpHelper::tFileHeader BmpFileHeader = {};
		BmpFileHeader.Offset = sizeof(MagicString) + sizeof(BmpHelper::tFileHeader) + sizeof(BmpHelper::tInfoHeader);
		ASSERT(BmpFileHeader.Offset == 54);
		BmpFileHeader.Size = BmpFileHeader.Offset + 12 + 4*SizeX*SizeY;
		if (!fwrite(&BmpFileHeader, sizeof(BmpFileHeader), 1, StreamOut)) {
			fprintf(gLogStream, "ERROR: fail writing to output!!!");
			exit(-1);
		}
		
		BmpHelper::tInfoHeader BmpInfoHeader = {};
		BmpInfoHeader.BitCount = 32;
		BmpInfoHeader.Compression = 3;
		BmpInfoHeader.SizeX = SizeX;
		BmpInfoHeader.SizeY = SizeY;
		BmpInfoHeader.HeaderSize = sizeof(BmpInfoHeader);
		ASSERT(BmpInfoHeader.HeaderSize == 40);
		BmpInfoHeader.SizeImage = 4*SizeX*SizeY;
		if (!fwrite(&BmpInfoHeader, sizeof(BmpInfoHeader), 1, StreamOut)) {
			fprintf(gLogStream, "ERROR: fail writing to output!!!");
			exit(-1);
		}
		
		BmpHelper::PutLayerToBmp32(StreamOut, BmpInfoHeader, &Layers[0][InitLevel], &Layers[1][InitLevel], &Layers[2][InitLevel], &Layers[3][InitLevel]);
	}
	
	STOP_CLOCK(DeCode);
}

#ifdef TEST
	static inline
	//=============================================================================
	void Test(
	//=============================================================================
	) {
		HaarWavelet::Test();
		HilbertCurve::Test();
		FibonacciCode::Test();
		Huffman::Test();
//TODO		Layer::Test();
	}
#endif


//=============================================================================
int main(
	int      ArgCount,
	tChar8** Args
//=============================================================================
) {
	START_CLOCK(main);
	
	gLogStream = stderr;
	
	FibonacciCode::Init(FibonacciCode::gFibArray, _countof(FibonacciCode::gFibArray));
	
#	ifdef TEST
		if (ArgCount > 1 && strncmp(Args[1], "-t", 2) == 0) {
			
			Test();
			if (ArgCount == 2) {
				exit(-1);
			}
			ArgCount -= 1;
			Args = &Args[1];
		}
#	endif
	
	if (ArgCount < 2) {
		printf("use -eN for encode or -d for decode");
		exit(-1);
	}
	
	auto Flag = Args[1];
	auto StreamIn = stdin; // default
	if (ArgCount > 2) {
		auto FileInName = Args[2];
		
		StreamIn = fopen(FileInName, "rb");
		if (!StreamIn) {
			fprintf(gLogStream, "ERROR: File '%s' not found!!!", FileInName);
			exit(-1);
		}
	}
	
	auto StreamOut = stdout; // default
	if (ArgCount > 3) {
		auto FileInName = Args[3];
		
		StreamOut = fopen(FileInName, "wb");
		if (!StreamIn) {
			fprintf(gLogStream, "ERROR: File '%s' not found!!!", FileInName);
			exit(-1);
		}
	}
	
	if (strncmp(Flag, "-e", 2) == 0) {
		BufferdStream::tWriter BufferdStreamOut;
		BufferdStream::Init(&BufferdStreamOut, StreamOut);
		EnCode(StreamIn, &BufferdStreamOut, atoi(&Flag[2]));
		BufferdStream::Close(&BufferdStreamOut);
	} else
	if (strncmp(Flag, "-d", 2) == 0) {
		BufferdStream::tReader BufferdStreamIn;
		BufferdStream::Init(&BufferdStreamIn, StreamIn);
		DeCode(&BufferdStreamIn, StreamOut);
		BufferdStream::Close(&BufferdStreamIn);
	}
	
	STOP_CLOCK(main);
#	ifdef VS
		PRINT_CLOCKS(fopen("C:\\Projekte\\ImgCodec2015\\bin\\timer.txt", "wb"));
#	else
		PRINT_CLOCKS(gLogStream);
#	endif
	return 0;
}
