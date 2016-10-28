#pragma warning(push, 0)

#include <stdio.h> // fopen(...), fclose(...), fread(...), fwrite(...), fflush(...), fprintf(...)
#include <string.h> // strncmp(...)
#include <stdlib.h> // exit(...)

#ifdef VS
#	include <malloc.h> // alloca(...)
#	define LL "ll"
#endif

#ifndef VS
#	define LL "I64"
#endif

#ifdef TEST
#	ifndef DEBUG
#		define DEBUG
#	endif
#endif

#pragma warning( pop )

// TODO-List:
// ==========
// - remove all arrays by tArray
// - StackBased-MemoryManager (replace malloc and free)
// - DependencyInjection for MemoryManagment
// - Huffman vs Fibonacci
// - Parallel
// - ContainerFormat

using tBool = bool;

using tNat1 = unsigned char;
using tNat8 = unsigned char;
using tNat16 = unsigned short;
using tNat32 = unsigned int;
using tNat64 = unsigned long long;

using tInt8 = signed char;
using tInt16 = signed short;
using tInt32 = signed int;
using tInt64 = signed long long;

using tReal32 = float;
using tReal64 = double;

const tBool cTrue = true;
const tBool cFalse = false;

using tChar8 = char;

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
#	ifndef TEST
#		define TEST
#	endif
#	ifndef TIMER
#		define TIMER
#	endif
	
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
#				ifndef VS
#					pragma GCC diagnostic push
#					pragma GCC diagnostic ignored "-Wdiv-by-zero"
#				endif
				LineNr /= 0;
				exit(-1);
#				ifndef VS
#					pragma GCC diagnostic pop
#				endif
			}
		}
		
	}
#else
#	define ASSERT(Cond)
#endif

//=============================================================================
template <typename g>
static inline
g Abs(
	g Value
//=============================================================================
) {
	return (Value < 0) ? -Value : Value;
}

//=============================================================================
static inline
tInt8 Sign(
	tInt64 Value
//=============================================================================
) {
	if (Value == 0) { return 0; }
	return Value >> 63;
}

//=============================================================================
template <typename g>
static inline
g Min(
	g A,
	g B
//=============================================================================
) {
	return (A < B) ? A : B;
}

//=============================================================================
template <typename g>
static inline
g Max(
	g A,
	g B
//=============================================================================
) {
	return (A > B) ? A : B;
}

//=============================================================================
static inline
tNat8 UsedBits(
	tNat64 Nat
//=============================================================================
) {
	auto Bits = (tNat8)0;
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
struct tArray final {
	g*     Values;
	tNat32 Size;
	tNat8  _[4];

	tArray(
	) {
		this->Size = 0;
		this->Values = NULL;
	}
	
	tArray(
		tNat32 Size,
		g* Values
	) {
		this->Size = Size;
		this->Values = Values;
	}
	
	g& operator[] (tNat64);
};

//=============================================================================
template <typename g>
inline
g& tArray<g>::operator[] (
	tNat64 Index
//=============================================================================
) {
	ASSERT(Index < this->Size);
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
tArray<g> MemAlloc(
	tArray<tNat8>& Buffer,
	tNat32         Size
//=============================================================================
) {
	auto Bytes = (tNat32)(sizeof(g) * Size);
	
	ASSERT(Bytes <= Buffer.Size);
	
	tArray<g> Result;
	Result.Size = Size;
	Result.Values = (Size == 0) ? NULL : (g*)Buffer.Values;
	
	Buffer.Values += Bytes;
	Buffer.Size -= Bytes;
	
	return Result;
}

//=============================================================================
template<typename g>
static inline
void MemZero(
	g* Buffer,
	tNat32 Count
//=============================================================================
) {
	// TODO: performance
	for (auto I = (tNat32)(sizeof(g) * Count); I --> 0; ) {
		((tNat8*)Buffer)[I] = 0;
	}
}

//=============================================================================
template<typename g>
static inline
void InPlaceQuickSort(
	tArray<g>            Nodes,
	tFunc2<g, g, tInt8>* Comp
//=============================================================================
) {
	if (Nodes.Size <= 1) {
		return;
	}
	
	auto Pivot = Nodes[Nodes.Size >> 1];
	auto Low   = (tNat32)0;
	auto Hight = (tNat32)(Nodes.Size - 1);
	
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
		auto Comp = (tFunc2<tNat32, tNat32, tInt8>*) [](tNat32 A, tNat32 B) { return Sign((tInt64)A - (tInt64)B); };
		while (CountIndex < Array.Size) {
			auto Count = Array[CountIndex];
			auto T = Take(Skip(Array, CountIndex + 1        ), Count);
			auto R = Take(Skip(Array, CountIndex + 1 + Count), Count);
			InPlaceQuickSort(T, Comp);
			for (auto I = Count; I --> 0; ) {
				ASSERT(T[I] == R[I]);
			}
			CountIndex += 2*Count + 1;
		}
	}
	
	//=============================================================================
	static inline
	void TestArray(
	//=============================================================================
	) {
		TestSwap();
		TestSort();
	}
#endif

#ifdef TIMER
#	define MEASURE_BLOCK(ClockId) Clock::tMeasure ClockId(Clock::Id::ClockId);
#	define START_CLOCK(ClockId) Clock::Start(Clock::Id::ClockId);
#	define STOP_CLOCK(ClockId) Clock::Stop(Clock::Id::ClockId);
#	define PRINT_CLOCKS(Stream) Clock::print((Stream));
	
	namespace Clock {
		
		using tClock = struct {
			tNat64 Begin    = 0;
			tNat64 Duration = 0;
			tNat64 Count    = 0;
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
				WriteLayer_ShowHistogram,
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
				
				MaxClocks
			};
		}
		
		tClock gClocks[Id::MaxClocks] = {};
		
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
			if (gClocks[Id].Begin != 0) {
				gClocks[Id].Duration += clock() - gClocks[Id].Begin;
				gClocks[Id].Begin = 0;
				gClocks[Id].Count += 1;
			}
		}
		
		class tMeasure final {
			tNat32 _Id = 0;
			
			public:
			//=============================================================================
			inline
			tMeasure(
				tNat32 aId
			//=============================================================================
			) {
				_Id = aId;
				Start(_Id);
			}
			
			//=============================================================================
			inline
			~tMeasure(
			//=============================================================================
			) {
				Stop(_Id);
			}
		};
		
		//=============================================================================
		static inline
		void print(
			tStream* Stream
		//=============================================================================
		) {
			static const char* Names[Id::MaxClocks] = {};
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
			ADD_CLOCK(WriteLayer_ShowHistogram);
			ADD_CLOCK(WriteByte_fwrite);
			ADD_CLOCK(WriteNat_);
			ADD_CLOCK(ReadByte_fread);
#			undef ADD_CLOCK
			
			const auto c1k = (tNat64)1000;
			const auto c10k = (tNat64)10000;
			
			for (auto I = (tNat32)0; I < Id::MaxClocks; I += 1) {
				auto Count = Clock::gClocks[I].Count;
				if (Count != 0) {
					auto Duration = Clock::gClocks[I].Duration;
					auto Mean = Duration / Count;
					
					auto CountUnit = "";
					if (Count > c10k) { Count /= c1k; CountUnit = "k"; }
					if (Count > c10k) { Count /= c1k; CountUnit = "M"; }
					
					auto DuratuonUnit = "";
					if (Duration > c10k) { Duration /= c1k; DuratuonUnit = "k"; }
					if (Duration > c10k) { Duration /= c1k; DuratuonUnit = "M"; }
					
					auto MeanUnit = "";
					if (Mean > c10k) { Mean /= c1k; MeanUnit = "k"; }
					if (Mean > c10k) { Mean /= c1k; MeanUnit = "M"; }
					
					if (Names[I]) {
						fprintf(Stream, "ClockId %s: ", Names[I]);
					} else {
						fprintf(Stream, "ClockId %d: ", I);
					}
					fprintf(Stream, "%" LL "u", Count);
					fprintf(Stream, "%s ", CountUnit);
					fprintf(Stream, "* %" LL "u", Mean);
					fprintf(Stream, "%s = ", MeanUnit);
					fprintf(Stream, "%" LL "u", Duration);
					fprintf(Stream, "%s\n", DuratuonUnit);
				}
			}
		}
		
	}
#else
#	define MEASURE_BLOCK(ClockId)
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
	Result.Size = (tNat32)Count;
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
//=============================================================================
template <typename g>
g* HeapAlloc(
//=============================================================================
) {
	return (g*)malloc(sizeof(g));
}

//=============================================================================
template <typename g>
void HeapFree(
	g* Memory
//=============================================================================
) {
	free(Memory);
}

namespace ByteStream {
	using tReader = struct {
		void*                 Env;
		tFunc1<void*, tNat8>* ReadByte;
		tFunc1<void*, void>*  Close;
	};
	
	//=============================================================================
	static inline
	tReader Reader(
		void*                 Env,
		tFunc1<void*, tNat8>* ReadByte,
		tFunc1<void*, void>*  Close
	//=============================================================================
	) {
		tReader Result;
		Result.Env = Env;
		Result.ReadByte = ReadByte;
		Result.Close = Close;
		return Result;
	}
	
	using tWriter = struct {
		void*                       Env;
		tFunc2<void*, tNat8, void>* WriteByte;
		tFunc1<void*, void>*        Close;
	};

	//=============================================================================
	static inline
	tWriter Writer(
		void*                       Env,
		tFunc2<void*, tNat8, void>* WriteByte,
		tFunc1<void*, void>*        Close
	//=============================================================================
	) {
		tWriter Result;
		Result.Env = Env;
		Result.WriteByte = WriteByte;
		Result.Close = Close;
		return Result;
	}
	
	//=============================================================================
	static inline
	tNat8 Read (
		ByteStream::tReader* Reader
	//=============================================================================
	) {
		return Reader->ReadByte(Reader->Env);
	}
	
	//=============================================================================
	static inline
	void Close (
		ByteStream::tReader* Reader
	//=============================================================================
	) {
		Reader->Close(Reader->Env);
	}
	
	//=============================================================================
	static inline
	void Write (
		ByteStream::tWriter* Writer,
		tNat8                Byte
	//=============================================================================
	) {
		Writer->WriteByte(Writer->Env, Byte);
	}
	
	//=============================================================================
	static inline
	void Close (
		ByteStream::tWriter* Writer
	//=============================================================================
	) {
		Writer->Close(Writer->Env);
	}
}

namespace MemoryStream {
	using tReader = struct {
		tArray<tNat8> Bytes;
		tNat32        Count;
		tNat8         _[4];
	};
	
	//=============================================================================
	static inline
	ByteStream::tReader Reader(
		tArray<tNat8> Bytes
	//=============================================================================
	) {
		auto Env = HeapAlloc<tReader>();
		Env->Bytes = Bytes;
		Env->Count = 0;
		return ByteStream::Reader(
			(void*)Env,
			[](void* Env) {
				auto Reader = (tReader*)Env;
				ASSERT(Reader->Count < Reader->Bytes.Size);
				auto Byte = Reader->Bytes[Reader->Count];
				Reader->Count += 1;
				return Byte;
			},
			[](void* Env) {
				HeapFree((tReader*)Env);
			}
		);
	}
	
	using tWriter = struct {
		tArray<tNat8> Bytes;
		tNat32        Count;
		tNat8         _[4];
	};
	
	//=============================================================================
	static inline
	ByteStream::tWriter Writer(
		tArray<tNat8> Bytes
	//=============================================================================
	) {
		auto Env = HeapAlloc<tReader>();
		Env->Bytes = Bytes;
		Env->Count = 0;
		return ByteStream::Writer(
			(void*)Env,
			[](void* Env, tNat8 Byte) {
				auto Writer = (tWriter*)Env;
				ASSERT(Writer->Count < Writer->Bytes.Size);
				Writer->Bytes[Writer->Count] = Byte;
				Writer->Count += 1;
			},
			[](void* Env) {
				HeapFree((tWriter*)Env);
			}
		);
	}
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
//		MEASURE_BLOCK(WriteByte);
		
		const auto cBufferSize = (tNat32)sizeof(Stream->Buffer);
		auto Pos = Stream->Pos;
		Stream->Buffer[Pos] = Byte;
		Pos = (Pos + 1) % cBufferSize;
		if (Pos == 0) {
//			MEASURE_BLOCK(WriteByte_fwrite);
			fwrite(Stream->Buffer, cBufferSize, 1, Stream->Stream);
		}
		Stream->Pos = Pos;
	}
	
	static inline
	//=============================================================================
	tNat8 ReadByte(
		tReader* Stream
	//=============================================================================
	) {
		MEASURE_BLOCK(ReadByte);
		
		const auto cBufferSize = (tNat32)sizeof(Stream->Buffer);
		auto Pos = Stream->Pos;
		if (Pos == 0) {
//			START_CLOCK(ReadByte_fread);
			auto Bytes = (tNat32)fread(Stream->Buffer, 1, cBufferSize, Stream->Stream);
//			STOP_CLOCK(ReadByte_fread);
			if (Bytes == 0) {
				return 0; // TODO: Fehlerbehandlung
			}
			if (Bytes != cBufferSize) {
				auto Des = cBufferSize;
				auto Src = Bytes;
				while (Src > 0) {
					Des -= 1;
					Src -= 1;
					Stream->Buffer[Des] = Stream->Buffer[Src];
				}
				Pos = cBufferSize - Bytes;
			}
		}
		Stream->Pos = (Pos + 1) % cBufferSize;
		auto Result =  Stream->Buffer[Pos];
		
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
//		START_CLOCK(WriteByte_fwrite);
		fwrite(Stream->Buffer, Stream->Pos, 1, Stream->Stream);
//		STOP_CLOCK(WriteByte_fwrite);
		Stream->Pos = 0;
		Stream->Stream = NULL;
	}
}

using iBitReader = struct {
	tFunc1<void*, tNat1>* ReadBit;
	tFunc1<void*, void>*  Close;
	void*                 Context;
};

//=============================================================================
static inline
tNat1 ReadBit(
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
	tFunc2<void*, tNat1, void>* WriteBit;
	tFunc1<void*, void>*        Close;
	void*                       Context;
};

//=============================================================================
static inline
void WriteBit(
	iBitWriter* BitWriter,
	tNat1       Bit
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
		ByteStream::tReader* Stream;
		tNat8                Bits = 0;
		tNat8                Byte = 0;
		tNat8                _[6];
	};
	
	using tWriter = struct {
		ByteStream::tWriter* Stream;
		tNat8                Bits = 0;
		tNat8                Byte = 0;
		tNat8                _[6];
	};
	
	//=============================================================================
	static inline
	void Init(
		tWriter*             BitStream,
		ByteStream::tWriter* Stream
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
		if (BitStream->Bits != 0) {
			ByteStream::Write(BitStream->Stream, BitStream->Byte);
		}
		*BitStream = {};
	}
	
	//=============================================================================
	static inline
	void WriteBit(
		tWriter* BitStream,
		tNat1    Bit
	//=============================================================================
	) {
//		MEASURE_BLOCK(WriteBit);
		
		ASSERT((Bit & ~1) == 0);
		
		auto Bits = BitStream->Bits;
		
		if (Bits == 0) {
			BitStream->Byte = Bit;
		} else {
			BitStream->Byte |= Bit << Bits;
		}
		
		Bits += 1;
		Bits &= 7;
		
		if (Bits == 0) {
			Write(BitStream->Stream, BitStream->Byte);
		}
		BitStream->Bits = Bits;
	}
	
	//=============================================================================
	static inline
	void Init(
		tReader*             BitStream,
		ByteStream::tReader* Stream
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
	tNat1 ReadBit(
		tReader* BitStream
	//=============================================================================
	) {
//		MEASURE_BLOCK(ReadBit);
		
		auto Bits = BitStream->Bits;
		auto Byte = (tNat8)0;
		if (Bits == 0) {
			Byte = ByteStream::Read(BitStream->Stream);
		} else {
			Byte = BitStream->Byte >> 1;
		}
		
		Bits = (Bits + 1) & 7;
		
		BitStream->Bits = Bits;
		BitStream->Byte = Byte;
		
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
		Interface.WriteBit = (tFunc2<void*, tNat1, void>*)WriteBit;
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
		Interface.ReadBit = (tFunc1<void*, tNat1>*)ReadBit;
		Interface.Close   = (tFunc1<void*, void>*)(tFunc1<tReader*, void>*)Close;
		return Interface;
	}
}

namespace ArithmeticBitStream {
	using tWriter = struct {
		iBitWriter* BitStream;
		tNat32 UnknowBits = 0;
		tNat32 Min        = 0;
		tNat32 Range      = 0;
		tNat32 Count0     = 0;
		tNat32 Count      = 0;
		tNat32 Count0_Old = 0;
	};
	
	using tReader = struct {
		iBitReader* BitStream;
		tNat32 UnknowBits = 0;
		tNat32 Value      = 0;
		tNat32 Min        = 0;
		tNat32 Range      = 0;
		tNat32 Count0     = 0;
		tNat32 Count      = 0;
		tNat32 Count0_Old = 0;
		tNat8  _[4]       = {0, 0, 0, 0}; // 64 bit alignment
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
		tNat1    Bit
	//=============================================================================
	) {
//		MEASURE_BLOCK(ArithmeticBitStream_WriteBit);
		
		//const tNat32 Full = 1 << 16;
		const auto cHalf = (tNat32)(1 << 15);
		const auto cQuad = (tNat32)(1 << 14);
		
		auto Delta = (Stream->Count0_Old * Stream->Range) >> 8;
		
		if (Bit == 0) {
			Stream->Range  = Delta;
		} else {
			Stream->Min   += Delta;
			Stream->Range -= Delta;
		}
		
		for(;;) {
			if ((Stream->Min + Stream->Range) - 1 < cHalf) {
				Stream->Range <<= 1;
				Stream->Min    -= 0;
				Stream->Min   <<= 1;
				
				WriteBit(Stream->BitStream, 0);
				while (Stream->UnknowBits != 0) {
					Stream->UnknowBits -= 1;
					WriteBit(Stream->BitStream, 1);
				}
			} else if (Stream->Min >= cHalf) {
				Stream->Range <<= 1;
				Stream->Min    -= cHalf;
				Stream->Min   <<= 1;
				
				WriteBit(Stream->BitStream, 1);
				while (Stream->UnknowBits != 0) {
					Stream->UnknowBits -= 1;
					WriteBit(Stream->BitStream, 0);
				}
			} else if (
				Stream->Min >= cQuad &&
				(Stream->Min + Stream->Range) - 1 < 3*cQuad
			) {
				Stream->Range <<= 1;
				Stream->Min    -= cQuad;
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
			Stream->Count0_Old = Stream->Count0;
		}
	}
	
	//=============================================================================
	static inline
	tNat1 ReadBit(
		tReader* Stream
	//=============================================================================
	) {
//		MEASURE_BLOCK(ArithmeticBitStream_ReadBit);
		
		//const tNat32 Full = 1 << 16;
		const auto cHalf = (tNat32)(1 << 15);
		const auto cQuad = (tNat32)(1 << 14);
		
		auto Bit = (tNat32)0;
		
		if (Stream->UnknowBits != 0) {
			Bit = Stream->Value;
			Stream->UnknowBits -= 1;
		} else {
			for (;;) {
				if ((Stream->Min + Stream->Range) - 1 < cHalf) {
					Bit = 0;
					Stream->Value = 1;
					break;
				} else if (Stream->Min >= cHalf) {
					Bit = 1;
					Stream->Value = 0;
					break;
				} else if (
					Stream->Min >= cQuad &&
					(Stream->Min + Stream->Range) - 1 < 3*cQuad
				) {
					Stream->UnknowBits += 1;
					Stream->Range <<= 1;
					Stream->Min    -= cQuad;
					Stream->Min   <<= 1;
				} else {
					if (ReadBit(Stream->BitStream) == 0) {
						Stream->Range <<= 1;
						Stream->Min    -= 0;
						Stream->Min   <<= 1;
					} else {
						Stream->Range <<= 1;
						Stream->Min    -= cHalf;
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
			Stream->Count0_Old = Stream->Count0;
			Stream->Count0 >>= 1;
		}
		
		auto Delta = (Stream->Count0_Old * Stream->Range) >> 8;
		
		if (Bit == 0) {
			Stream->Range  = Delta;
		} else {
			Stream->Min   += Delta;
			Stream->Range -= Delta;
		}
		
		return (tNat1)Bit;
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
		Interface.ReadBit = (tFunc1<void*, tNat1>*)(tFunc1<tReader*, tNat1>*)ReadBit;
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
		Interface.WriteBit = (tFunc2<void*, tNat1, void>*)(tFunc2<tWriter*, tNat1, void>*)WriteBit;
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
//	MEASURE_BLOCK(WriteNat_);
	
	for (auto I = (tNat8)16; I --> 0; ) {
		auto Bit = (tNat1)((Value >> I) & 1);
		WriteBit(BitWriter, Bit);
	}
}

//=============================================================================
static inline
tNat64 ReadNat_(
	iBitReader* BitReader
//=============================================================================
) {
//	MEASURE_BLOCK(ReadNat_);
	
	auto Result = (tNat64)0;
	for (auto I = (tNat8)16; I --> 0; ) {
		Result <<= 1;
		Result |= ReadBit(BitReader);
	}
	
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
//		MEASURE_BLOCK(Elias_Write);
		
		Value += 1;
		
		auto DataBits = (tNat8)1;
		for (auto Temp = Value >> 1; Temp > 0; Temp >>= 1) {
			WriteBit(BitWriter, 0);
			DataBits += 1;
		}
		
		while (DataBits --> 0) {
			WriteBit(BitWriter, (Value >> DataBits) & 1);
		}
	}
	
	//=============================================================================
	static inline
	tInt64 Read(
		iBitReader* BitReader
	//=============================================================================
	) {
//		MEASURE_BLOCK(Elias_Read);
		
		auto DataBits = (tNat8)1;
		auto Bit = ReadBit(BitReader);
		while (Bit == 0) {
			DataBits += 1;
			Bit = ReadBit(BitReader);
		}
		
		auto Value = (tInt64)1;
		while (DataBits --> 0) {
			Value <<= 1;
			Value  |= ReadBit(BitReader);
		}
		
		return Value - 1;
	}
}

namespace FibonacciCode {
	static tNat64 gFibArray[50];
	
	//=============================================================================
	static inline
	void Init(
		tNat64* FibArray,
		tNat8  MaxFibBits
	//=============================================================================
	) {
		auto F1 = (tNat64)1;
		auto F2 = (tNat64)1;

		ASSERT(MaxFibBits == 50)
		
		for (auto FibBits = (tNat8)0; FibBits <= MaxFibBits; FibBits += 1) {
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
//		MEASURE_BLOCK(NatToFib);
		
		auto FibBits = UsedBits(Nat);
		FibBits += FibBits >> 1;
		
		auto Fib = (tNat64)0;
		while (FibBits --> 0) {
			auto FibI = FibArray[FibBits];
			auto Bit = (tNat1)((Nat >= FibI) ? 1 : 0);
			if (Bit == 1) {
				Nat -= FibI;
			}
			Fib <<= 1;
			Fib |= Bit;
		}
		
		return Fib;
	}
	
	//=============================================================================
	static inline
	tNat64 FibToNat(
		tNat64 Fib
	//=============================================================================
	) {
//		MEASURE_BLOCK(FibToNat);
		
		auto Nat = (tNat64)0;
		
		auto F1 = (tNat64)0;
		auto F2 = (tNat64)1;
		for (auto I = (tNat8)0; Fib != 0; I += 1, Fib >>= 1) {
			auto F = F1 + F2;
			
			F1 = F2;
			F2 = F;
			
			if ((Fib & 1) != 0) {
				Nat += F;
			}
		}
		
		return Nat;
	}
	
	//=============================================================================
	static inline
	void Write(
		iBitWriter* BitWriter,
		tNat64      Value
	//=============================================================================
	) {
//		MEASURE_BLOCK(Fib_Write);
		auto Fib = NatToFib(gFibArray, Value + 1);
		
		while (Fib != 0) {
			WriteBit(BitWriter, Fib & 1);
			Fib >>= 1;
		}
		
		WriteBit(BitWriter, 1);
	}
	
	//=============================================================================
	static inline
	tNat64 Read(
		iBitReader* BitReader
	//=============================================================================
	) {
//		MEASURE_BLOCK(Fib_Read);
		
		auto Fib = (tNat64)ReadBit(BitReader);
		for (auto Pos = (tNat8)1; cTrue; Pos += 1) {
			auto Bit = (tNat64)ReadBit(BitReader);
			if ((Bit & (Fib >> (Pos - 1))) != 0) { break; }
			Fib |= Bit << Pos;
		}
		
		auto Result = FibToNat(Fib) - 1;
		
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

namespace Huffman {
	//=============================================================================
	static inline
	void CalculateBitCount(
		tArray<tNat32>& SortedHistogramm,
		tArray<tNat8>& BitCounts         // OUT
	//=============================================================================
	) {
		ASSERT(SortedHistogramm.Size <= (1 << 17))
		
		// init
		auto HistogrammSize = SortedHistogramm.Size;
		auto NodeCount = 2*HistogrammSize - 1;
		
		auto NodeIds     = HeapAlloc<tNat32>(NodeCount);
		auto ParentNodes = HeapAlloc<tNat32>(NodeCount);
		auto Counts      = HeapAlloc<tNat32>(NodeCount);
		
		for (auto NodeId = NodeCount; NodeId --> 0; ) {
			NodeIds[NodeId] = NodeId;
			ParentNodes[NodeId] = 0;
			Counts[NodeId] = (NodeId < HistogrammSize) ? SortedHistogramm[NodeId] : 0;
		}
		
		{ // build tree
			auto I = (tNat32)0;
			auto NewNodeId = HistogrammSize;
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
		for (auto I = NodeCount - 1; I --> 0; ) {
			auto NodeId       = NodeIds[I];
			auto ParentNodeId = ParentNodes[NodeId];
			Counts[NodeId]    = Counts[ParentNodeId] + 1;
			if (NodeId < HistogrammSize) {
				BitCounts[NodeId] = (tNat8)Counts[NodeId];
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
		tArray<tNat8>&  BitCounts,
		tArray<tNat32>& BitsArray  // OUT
	//=============================================================================
	) {
		auto Bits = (tNat32)0;
		auto BitCount = (tNat8)0;
		for (auto NodeId = BitCounts.Size; NodeId --> 0; ) {
			auto NewBitCount = BitCounts[NodeId];
			Bits <<= NewBitCount - BitCount;
			BitCount = NewBitCount;
			
			BitsArray[NodeId] = Bits;
			
			Bits += 1;
		}
	}
	
	using tCount_NodeIndex = struct {
		tNat32 Count;
		tNat32 NodeIndex;
	};
	
	//=============================================================================
	static inline
	tInt8 CompareCounts(
		tCount_NodeIndex A,
		tCount_NodeIndex B
	//=============================================================================
	) {
		return Sign((tInt64)A.Count - (tInt64)B.Count);
	}
	
	//=============================================================================
	static inline
	void Init(
		tArray<tNat32>& Counts,
		tArray<tNat8>&  BitCounts, // Out
		tArray<tNat32>& BitsArray  // OUT
	//=============================================================================
	) {
		auto Count = Counts.Size;
		auto SortedCounts = HeapAlloc<tNat32>(Count);
		auto SortedIndex = HeapAlloc<tNat32>(Count);
		{
			auto ForSortArray = HeapAlloc<tCount_NodeIndex>(Count);
			for (auto I = Count; I --> 0; ) {
				ForSortArray[Count].Count = Counts[Count];
				ForSortArray[Count].NodeIndex = Count;
			}
			InPlaceQuickSort(ForSortArray, CompareCounts);
			for (auto I = Count; I --> 0; ) {
				SortedCounts[I] = ForSortArray[Count].Count;
				SortedIndex[I] = ForSortArray[Count].NodeIndex;
			}
			HeapFree(ForSortArray);
		}
		{
			auto SortedBitCounts = HeapAlloc<tNat8>(Count);
			auto SortedBitsArray = HeapAlloc<tNat32>(Count);
			CalculateBitCount(SortedCounts, SortedBitCounts);
			CalculateBits(SortedBitCounts, SortedBitsArray);
			for (auto I = Count; I --> 0; ) {
				auto J = SortedIndex[I];
				BitCounts[J] = SortedBitCounts[I];
				BitsArray[J] = SortedBitsArray[I];
			}
			HeapFree(SortedBitsArray);
			HeapFree(SortedBitCounts);
		}
		HeapFree(SortedCounts);
		HeapFree(SortedIndex);
	}
	
	struct tWriter final {
		iBitWriter*    BitWriter;
		tArray<tNat8>  BitCounts;
		tArray<tNat32> Bits;
	};
	
	//=============================================================================
	static inline
	tWriter Writer(
		iBitWriter*    BitWriter,
		tArray<tNat8>  BitCounts,
		tArray<tNat32> Bits
	//=============================================================================
	) {
		tWriter Result;
		Result.BitWriter = BitWriter;
		Result.BitCounts = BitCounts;
		Result.Bits = Bits;
		return Result;
	}
	
	//=============================================================================
	static inline
	void Write(
		tWriter* Writer,
		tNat32   Id
	//=============================================================================
	) {
		auto BitCount = Writer->BitCounts[Id];
		auto Bits = Writer->Bits[Id];
		while (BitCount --> 0) {
			WriteBit(Writer->BitWriter, Bits & 1);
			Bits >>= 1;
		}
	}
	
	struct tReader final {
		iBitReader*    BitReader;
		tArray<tNat8>  BitCounts;
		tArray<tNat32> Bits;
	};
	
	//=============================================================================
	static inline
	tReader Read(
		iBitReader*    BitReader,
		tArray<tNat8>  BitCounts,
		tArray<tNat32> Bits
	//=============================================================================
	) {
		tReader Result;
		Result.BitReader = BitReader;
		Result.BitCounts = BitCounts;
		Result.Bits = Bits;
		return Result;
	}
	
	//=============================================================================
	static inline
	tNat32 Read(
		tReader* Reader
	//=============================================================================
	) {
		auto BitCount = (tNat8)0;
		auto Bits = (tNat32)0;
		for (auto I = Reader->BitCounts.Size; I --> 0; ) {
			if (Reader->BitCounts[I] > BitCount) {
				Bits <<= 1;
				Bits |= ReadBit(Reader->BitReader);
				BitCount += 1;
			}
			
			if (Bits == Reader->Bits[I]) {
				return I;
			}
		}
		ASSERT(cFalse);
	}
	
	//=============================================================================
	static inline
	void WriteHeader(
		iBitWriter*    BitWriter,
		tArray<tNat8>  BitCounts,
		tArray<tNat32> Bits
	//=============================================================================
	) {
		auto Count = BitCounts.Size;
		FibonacciCode::Write(BitWriter, Count);
		auto LastBitCount = (tNat8)0;
		for (auto I = (tNat32)0; I < Count; I += 1) {
			auto BitCount = BitCounts[I];
			auto Delta = (tInt16)BitCount - (tInt16)LastBitCount;
			auto IsNeg = (Delta >> 8) & 1;
			FibonacciCode::Write(BitWriter, (Delta << 1) ^ (0 - IsNeg));
			LastBitCount = BitCount;
			
			auto Bits_ = Bits[I];
			for (auto J = (tNat8)0; J < BitCount; J += 1) {
				WriteBit(BitWriter, (Bits_ >> (BitCount - J - 1)) & 1);
			}
		}
	}
	
	//=============================================================================
	static inline
	void ReadHeader(
		iBitReader*    BitReader,
		tArray<tNat8>  BitCounts, // OUT
		tArray<tNat32> Bits       // OUT
	//=============================================================================
	) {
		auto BitCount = (tNat8)0;
		auto Count = Bits.Size;
		for (auto I = (tNat32)0; I < Count; I += 1) {
			auto Value = FibonacciCode::Read(BitReader);
			auto IsNeg = Value & 1;
			BitCount += (tNat8)((Value >> 1) ^ (IsNeg ? 0xFF : 0));
			BitCounts[I] = BitCount;
			
			auto Bits_ = (tNat32)0;
			for (auto J = BitCount; J --> 0; ) {
				Bits_ <<= 1;
				Bits_ |= ReadBit(BitReader);
			}
			Bits[I] = Bits_;
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
			tNat8 BitCountsRef[] = { 6, 6, 5, 4, 3,  2,  2,  2 };
			tNat8 BitCountsRes[] = { 0, 0, 0, 0, 0,  0,  0,  0 };
			tNat32 BitsRes[] = {        0,        0,       0,      0,     0,    0,    0,    0 };
			tNat32 BitsRef[] = { 0b111111, 0b111110, 0b11110, 0b1110, 0b110, 0b10, 0b01, 0b00 };
			
			tArray<tNat32> InHistogramm = AsArray(Histogram);
			
			{
				tArray<tNat8> OutCounts = AsArray(BitCountsRes);
				CalculateBitCount(InHistogramm, OutCounts);
				for (auto I = _countof(Histogram); I --> 0; ) {
					ASSERT(BitCountsRes[I] == BitCountsRef[I]);
				}
				
				tArray<tNat32> OutBits = AsArray(BitsRes);
				CalculateBits(OutCounts, OutBits);
				for (auto I = _countof(Histogram); I --> 0; ) {
					ASSERT(BitsRes[I] == BitsRef[I]);
				}
			}
			
			auto Bytes = HeapAlloc<tNat8>(1<<10);
			{
				auto MemoryWriter = MemoryStream::Writer(Bytes);
				BitStream::tWriter BitWriter;
				BitStream::Init(&BitWriter, &MemoryWriter);
				auto BitWriter_ = BitStream::GetInterface(&BitWriter);
				Huffman::WriteHeader(&BitWriter_, AsArray(BitCountsRef), AsArray(BitsRef));
				FibonacciCode::Write(&BitWriter_, 7);
				auto HuffmanWriter = Huffman::Writer(&BitWriter_, AsArray(BitCountsRef), AsArray(BitsRef));
				for (auto I = 0; I < 7; I +=1) {
					Huffman::Write(&HuffmanWriter, I);
				}
				Close(&BitWriter);
				ByteStream::Close(&MemoryWriter);
			}
			for (auto I = _countof(BitCountsRef); I --> 0; ) {
				BitCountsRes[I] = 0;
				BitsRes[I] = 0;
			}
			{
				auto MemoryReader = MemoryStream::Reader(Bytes);
				BitStream::tReader BitReader;
				BitStream::Init(&BitReader, &MemoryReader);
				
				auto BitReader_ = BitStream::GetInterface(&BitReader);
				
				auto Count = FibonacciCode::Read(&BitReader_);
				auto OutCounts = HeapAlloc<tNat8>(Count);
				auto OutBits = HeapAlloc<tNat32>(Count);
				
				Huffman::ReadHeader(&BitReader_, OutCounts, OutBits);
				
				for (auto I = (tNat32)_countof(BitCountsRef); I --> 0; ) {
					ASSERT(BitCountsRef[I] == OutCounts[I]);
					ASSERT(BitsRef[I] == OutBits[I]);
				}
				
				HeapFree(OutBits);
				HeapFree(OutCounts);
				
				Close(&BitReader);
				
				ByteStream::Close(&MemoryReader);
			}
			HeapFree(Bytes);
		}
		
		//=============================================================================
		static inline
		void Test(
		//=============================================================================
		) {
			TestCalcBitCount();
		}
#	endif
}

namespace HaarWavelet {
	using tPair = struct {
		tInt16 A = 0;
		tInt16 B = 0;
	};
	
	//=============================================================================
	static inline
	tPair EnCode(
		tPair Src
	//=============================================================================
	) {
		tPair Res;
		Res.A = (Src.A + Src.B + 1) >> 1;
		Res.B = Src.A - Src.B;
		return Res;
	}
	
	//=============================================================================
	static inline
	tPair DeCode(
		tPair Src
	//=============================================================================
	) {
		tPair Res;
		Res.A = ((Src.A << 1) + Src.B) >> 1;
		Res.B = ((Src.A << 1) - Src.B) >> 1;
		return Res;
	}
	
	struct tQuad {
		tInt16 A = 0;
		tInt16 B = 0;
		tInt16 C = 0;
		tInt16 D = 0;
	};
	
	//=============================================================================
	static inline
	tQuad EnCode(
		tQuad Src
	//=============================================================================
	) {
		tQuad Res;
#		if 1
			Res.A = (Src.A + Src.B + Src.C + Src.D + 2) >> 2;
			Res.B =  Src.A - Src.B + Src.C - Src.D;
			Res.C =  Src.A + Src.B - Src.C - Src.D;
			Res.D =  Src.A - Src.B - Src.C + Src.D;
#		elif 0
			Res.A = (Src.A + Src.B + Src.C + Src.D + 1) >> 1;
			Res.B = (Src.A - Src.B + Src.C - Src.D + 0) >> 1;
			Res.C = (Src.A + Src.B - Src.C - Src.D + 0) >> 1;
			Res.D = (Src.A - Src.B - Src.C + Src.D + 0) >> 1;
			
			Res.D <<= 1;
			Res.D |= Res.A & 1;
			Res.A >>= 1;
#		elif 1
			Res.A = (Src.A + Src.B + Src.C + Src.D + 1) >> 1;
			Res.B = (Src.A - Src.B + Src.C - Src.D + 0) >> 1;
			Res.C = (Src.A + Src.B - Src.C - Src.D + 0) >> 1;
			Res.D = (Src.A - Src.B - Src.C + Src.D + 0) >> 1;
#		else
			tPair AB = { Src.A, Src.B };
			tPair CD = { Src.C, Src.D };
			
			auto AB_ = EnCode(AB);
			auto CD_ = EnCode(CD);
			
			tPair AC = { AB_.A, CD_.A };
			tPair BD = { AB_.B, CD_.B };
			
			auto AC_ = EnCode(AC);
			auto BD_ = EnCode(BD);
			
			Res.A = AC_.A;
			Res.B = BD_.A;
			Res.C = AC_.B;
			Res.D = BD_.B;
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
			Res.A = ((Src.A << 2) + 1 + Src.B + Src.C + Src.D) >> 2;
			Res.B = ((Src.A << 2) + 1 - Src.B + Src.C - Src.D) >> 2;
			Res.C = ((Src.A << 2) + 1 + Src.B - Src.C - Src.D) >> 2;
			Res.D = ((Src.A << 2) + 1 - Src.B - Src.C + Src.D) >> 2;
#		elif 0
			Src.A <<= 1;
			Src.A |= Src.D & 1;
			Res.D >>= 1;
			
			Res.A = (Src.A + Src.B + Src.C + Src.D + 1) >> 1;
			Res.B = (Src.A - Src.B + Src.C - Src.D + 0) >> 1;
			Res.C = (Src.A + Src.B - Src.C - Src.D + 0) >> 1;
			Res.D = (Src.A - Src.B - Src.C + Src.D + 0) >> 1;
#		elif 1
			Res.A = (Src.A + Src.B + Src.C + Src.D + 1) >> 1;
			Res.B = (Src.A - Src.B + Src.C - Src.D + 0) >> 1;
			Res.C = (Src.A + Src.B - Src.C - Src.D + 0) >> 1;
			Res.D = (Src.A - Src.B - Src.C + Src.D + 0) >> 1;
#		else
			tPair AC = { Src.A, Src.C };
			tPair BD = { Src.B, Src.D };
			
			auto AC_ = DeCode(AC);
			auto BD_ = DeCode(BD);
			
			tPair AB = { AC_.A, BD_.A };
			tPair CD = { AC_.B, BD_.B };
			
			auto AB_ = DeCode(AB);
			auto CD_ = DeCode(CD);
			
			Res.A = AB_.A;
			Res.B = AB_.B;
			Res.C = CD_.A;
			Res.D = CD_.B;
#		endif
		return Res;
	}
	
#	ifdef TEST
		//=============================================================================
		static inline
		void Test(
		//=============================================================================
		) {
			// 1-D
			tPair Src1D;
			
			// #1: X != Y -> EnCode(X) != EnCode(Y)
			for (Src1D.A = -2; Src1D.A < 2; Src1D.A += 1) {
				for (Src1D.B = -2; Src1D.B < 2; Src1D.B += 1) {
					auto Temp = EnCode(Src1D);
					auto Res  = DeCode(Temp);
					ASSERT(Src1D.A == Res.A);
					ASSERT(Src1D.B == Res.B);
				}
			}
			
			// #2: X != Y -> DeCode(X) != DeCode(Y)
			for (Src1D.A = -2; Src1D.A < 2; Src1D.A += 1) {
				for (Src1D.B = -2; Src1D.B < 2; Src1D.B += 1) {
					auto Temp = DeCode(Src1D);
					auto Res  = EnCode(Temp);
					ASSERT(Src1D.A == Res.A);
					ASSERT(Src1D.B == Res.B);
				}
			}
			
			// #1 & #2 -> EnCode & DeCode sind eineindeutige Abbildungen -> kompakte Kodierung
			
			// 2-D
			tQuad Src;
			
			// #1: X != Y -> EnCode(X) != EnCode(Y)
			for (Src.A = -4; Src.A < 4; Src.A += 1) {
				for (Src.B = -4; Src.B < 4; Src.B += 1) {
					for (Src.C = -4; Src.C < 4; Src.C += 1) {
						for (Src.D = -4; Src.D < 4; Src.D += 1) {
							auto Temp = EnCode(Src);
							auto Res  = DeCode(Temp);
							ASSERT(Res.A == Src.A);
							ASSERT(Res.B == Src.B);
							ASSERT(Res.C == Src.C);
							ASSERT(Res.D == Src.D);
						}
					}
				}
			}
			
			#if 0
			// #2: X != Y -> DeCode(X) != DeCode(Y)
			for (Src.A = -4; Src.A < 4; Src.A += 1) {
				for (Src.B = -4; Src.B < 4; Src.B += 1) {
					for (Src.C = -4; Src.C < 4; Src.C += 1) {
						for (Src.D = -4; Src.D < 4; Src.D += 1) {
							auto Temp = DeCode(Src);
							auto Res  = EnCode(Temp);
							ASSERT(Res.A == Src.A);
							ASSERT(Res.B == Src.B);
							ASSERT(Res.C == Src.C);
							ASSERT(Res.D == Src.D);
							auto X = 0;
						}
					}
				}
			}
			#endif
			// #1 -> Verlustfrei
			// #1 & #2 -> EnCode & DeCode sind eineindeutige Abbildungen -> kompakte Kodierung
		}
#	endif
}

namespace Layer {
	using tLevel = struct {
		tNat16         SizeXLeft;
		tNat16         SizeXRight;      //    Left  Right
		tNat16         SizeYTop;        //   +-----+---+
		tNat16         SizeYBottom;     //   |  S  | H | Top
		tArray<tInt16> S; // Summe      //   |     |   |
		tArray<tInt16> H; // Horizontal //   +-----+---+
		tArray<tInt16> V; // Vertical   //   |  V  | D | Bottom
		tArray<tInt16> D; // Diagonal   //   +-----+---+
		tNat16         Pitch;           // Pitch >= SizeXLeft >= SizeXRight; SizeYTop >= SizeYBottom
		tNat8          _[6];            // 64 bit alignment
	};
	
	//=============================================================================
	static
	tNat8 Init(
		tLevel*        Layers,
		tNat16         SizeX,
		tNat16         SizeY,
		tArray<tNat8>& Mem
	//=============================================================================
	) {
		auto InitLevel = (tNat8)0;
		{ // Berechne InitLevel + Initalisierung vom InitLevel
			auto X = SizeX;
			auto Y = SizeY;
			while (X > 1 || Y > 1) {
				X = (X + 1) >> 1;
				Y = (Y + 1) >> 1;
				InitLevel += 1;
			}
			auto Layer = &Layers[InitLevel];
			
			Layer->Pitch = (SizeX + 7) & ~7;
			Layer->SizeXLeft = SizeX;
			Layer->SizeXRight = 0;
			Layer->SizeYTop = SizeY;
			Layer->SizeYBottom = 0;
			
			Layer->S = MemAlloc<tInt16>(Mem, Layer->Pitch * Layer->SizeYTop);
			Layer->H = MemAlloc<tInt16>(Mem, 0);
			Layer->V = MemAlloc<tInt16>(Mem, 0);
			Layer->D = MemAlloc<tInt16>(Mem, 0);
		}
		
		{ // Initalisierung restlicher Levels
			auto X = (tNat16)SizeX;
			auto Y = (tNat16)SizeY;
			for (auto Level = InitLevel; Level --> 0; ) {
				auto Layer = &Layers[Level];
				
				Layer->SizeXRight  = (X + 0) >> 1;
				Layer->SizeYBottom = (Y + 0) >> 1;
				Layer->SizeXLeft   = (X + 1) >> 1;
				Layer->SizeYTop    = (Y + 1) >> 1;
				
				Layer->Pitch = (Layer->SizeXLeft + 7) & ~7;
				
				Layer->S = MemAlloc<tInt16>(Mem, Layer->Pitch * Layer->SizeYTop);
				Layer->H = MemAlloc<tInt16>(Mem, Layer->Pitch * Layer->SizeYTop);
				Layer->V = MemAlloc<tInt16>(Mem, Layer->Pitch * Layer->SizeYBottom);
				Layer->D = MemAlloc<tInt16>(Mem, Layer->Pitch * Layer->SizeYBottom);
				
				X = Layer->SizeXLeft;
				Y = Layer->SizeYTop;
			}
		}
		return InitLevel;
	}
	
	//=============================================================================
	static inline
	tInt16 DeltaQubicWavelet(
		tInt16 Last2,
		tInt16 Last1,
		tInt16 Curr,
		tInt16 Next1,
		tInt16 Next2
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
	tInt16 DeltaQubicWaveletPre2(
		tInt16 Curr,
		tInt16 Next1,
		tInt16 Next2
	//=============================================================================
	) {
		auto Last2 = (tInt16)(2*Curr - Next2);
		auto Last1 = (tInt16)(2*Curr - Next1);
		return DeltaQubicWavelet(Last2, Last1, Curr, Next1, Next2);
	}
	
	//=============================================================================
	static inline
	tInt16 DeltaQubicWaveletPre1(
		tInt16 Last1,
		tInt16 Curr,
		tInt16 Next1,
		tInt16 Next2
	//=============================================================================
	) {
		auto Last2 = (tInt16)(2*Last1 - Curr);
		return DeltaQubicWavelet(Last2, Last1, Curr, Next1, Next2);
	}
	
	//=============================================================================
	static inline
	tInt16 DeltaQubicWaveletPost1(
		tInt16 Last2,
		tInt16 Last1,
		tInt16 Curr,
		tInt16 Next1
	//=============================================================================
	) {
		auto Next2 = (tInt16)(2*Next1 - Curr);
		return DeltaQubicWavelet(Last2, Last1, Curr, Next1, Next2);
	}
	
	//=============================================================================
	static inline
	tInt16 DeltaQubicWaveletPost2(
		tInt16 Last2,
		tInt16 Last1,
		tInt16 Curr
	//=============================================================================
	) {
		auto Next2 = (tInt16)(2*Curr - Last2);
		auto Next1 = (tInt16)(2*Curr - Last1);
		return DeltaQubicWavelet(Last2, Last1, Curr, Next1, Next2);
	}
	
	const auto cPredDist = (tNat8)2;
	const auto cDFactor = (tNat8)2; // shift 2
	
	//=============================================================================
	static
	void SplitWithQubPred(
		tLevel* SrcLevel,
		tLevel* DesLevel
	//=============================================================================
	) {
		MEASURE_BLOCK(SplitWithQubPred);
		
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
		
		auto Last2_S = (tInt16)0;
		auto Last2_V = (tInt16)0;
		auto Last1_S = (tInt16)0;
		auto Last1_V = (tInt16)0;
		auto Curr0_S = (tInt16)0;
		auto Curr0_V = (tInt16)0;
		auto Next1_S = (tInt16)0;
		auto Next1_V = (tInt16)0;
		auto Next2_S = (tInt16)0;
		auto Next2_V = (tInt16)0;
		
		for (auto Y = -cPredDist; Y < SizeY0; Y++) {
			auto Y2 = Y << 1;
			for (auto X = -cPredDist; X < SizeX0; X++) {
				Last2_S = Last1_S;
				Last2_V = Last1_V;
				Last1_S = Curr0_S;
				Last1_V = Curr0_V;
				Curr0_S = Next1_S;
				Curr0_V = Next1_V;
				Next1_S = Next2_S;
				Next1_V = Next2_V;
				if (X + cPredDist < SizeX0) {
					if (Y + cPredDist < SizeY0) {
						auto X2 = X << 1;
						HaarWavelet::tQuad Arg;
						HaarWavelet::tQuad Res;
						if (Y2 + 2*cPredDist == SizeY - 1) {
							if (X2 + 2*cPredDist  == SizeX - 1) {
								Arg.A = Src[(tNat32)((X2 + 2*cPredDist + 0) + (Y2 + 2*cPredDist + 0) * SrcPitch)];
								Arg.B = Arg.A;
								Arg.C = Arg.A;
								Arg.D = Arg.A;
								
								Res = HaarWavelet::EnCode(Arg);
								
								DesS[(X + cPredDist) + (Y + cPredDist) * DesPitch] = Res.A;
							} else {
								Arg.A = Src[(tNat32)((X2 + 2*cPredDist + 0) + (Y2 + 2*cPredDist + 0) * SrcPitch)];
								Arg.B = Src[(tNat32)((X2 + 2*cPredDist + 1) + (Y2 + 2*cPredDist + 0) * SrcPitch)];
								Arg.C = Arg.A;
								Arg.D = Arg.B;
								
								Res = HaarWavelet::EnCode(Arg);
								
								DesS[(X + cPredDist) + (Y + cPredDist) * DesPitch] = Res.A;
								DesH[(X + cPredDist) + (Y + cPredDist) * DesPitch] = Res.B;
							}
						} else {
							if (X2 + 2*cPredDist == SizeX - 1) {
								Arg.A = Src[(tNat32)((X2 + 2*cPredDist + 0) + (Y2 + 2*cPredDist + 0) * SrcPitch)];
								Arg.B = Arg.A;
								Arg.C = Src[(tNat32)((X2 + 2*cPredDist + 0) + (Y2 + 2*cPredDist + 1) * SrcPitch)];
								Arg.D = Arg.C;
								
								Res = HaarWavelet::EnCode(Arg);
								
								DesS[(X + cPredDist) + (Y + cPredDist) * DesPitch] = Res.A;
								DesV[(X + cPredDist) + (Y + cPredDist) * DesPitch] = Res.C;
							} else {
								Arg.A = Src[(tNat32)((X2 + 2*cPredDist + 0) + (Y2 + 2*cPredDist + 0) * SrcPitch)];
								Arg.B = Src[(tNat32)((X2 + 2*cPredDist + 1) + (Y2 + 2*cPredDist + 0) * SrcPitch)];
								Arg.C = Src[(tNat32)((X2 + 2*cPredDist + 0) + (Y2 + 2*cPredDist + 1) * SrcPitch)];
								Arg.D = Src[(tNat32)((X2 + 2*cPredDist + 1) + (Y2 + 2*cPredDist + 1) * SrcPitch)];
								
								Res = HaarWavelet::EnCode(Arg);
								
								DesS[(X + cPredDist) + (Y + cPredDist) * DesPitch] = Res.A;
								DesH[(X + cPredDist) + (Y + cPredDist) * DesPitch] = Res.B;
								DesV[(X + cPredDist) + (Y + cPredDist) * DesPitch] = Res.C;
								DesD[(X + cPredDist) + (Y + cPredDist) * DesPitch] = Res.D;
							}
						}
					}
					
					if (Y < 0 || (SizeX & SizeY) < 8) {
						continue;
					}
					
					Next2_S = DesS[(X + cPredDist) + (Y + 0) * DesPitch];
					Next2_V = DesV[(X + cPredDist) + (Y + 0) * DesPitch];
					if (Y < SizeY1) {
						auto Temp = 0;
						if (Y == 0) {
							auto Curr0L_S = DesS[(X + cPredDist) + (Y + 0) * DesPitch];
							auto Next1L_S = DesS[(X + cPredDist) + (Y + 1) * DesPitch];
							auto Next2L_S = DesS[(X + cPredDist) + (Y + 2) * DesPitch];
							Temp = DeltaQubicWaveletPre2(Curr0L_S, Next1L_S, Next2L_S);
						} else
						if (Y == 1) {
							auto Last1L_S = DesS[(X + cPredDist) + (Y - 1) * DesPitch];
							auto Curr0L_S = DesS[(X + cPredDist) + (Y + 0) * DesPitch];
							auto Next1L_S = DesS[(X + cPredDist) + (Y + 1) * DesPitch];
							auto Next2L_S = DesS[(X + cPredDist) + (Y + 2) * DesPitch];
							Temp = DeltaQubicWaveletPre1(Last1L_S, Curr0L_S, Next1L_S, Next2L_S);
						} else
						if (Y + 2 < SizeY0) {
							auto Last2L_S = DesS[(X + cPredDist) + (Y - 2) * DesPitch];
							auto Last1L_S = DesS[(X + cPredDist) + (Y - 1) * DesPitch];
							auto Curr0L_S = DesS[(X + cPredDist) + (Y + 0) * DesPitch];
							auto Next1L_S = DesS[(X + cPredDist) + (Y + 1) * DesPitch];
							auto Next2L_S = DesS[(X + cPredDist) + (Y + 2) * DesPitch];
							Temp = DeltaQubicWavelet(Last2L_S, Last1L_S, Curr0L_S, Next1L_S, Next2L_S);
						} else
						if (Y + 2 == SizeY0) {
							auto Last2L_S = DesS[(X + cPredDist) + (Y - 2) * DesPitch];
							auto Last1L_S = DesS[(X + cPredDist) + (Y - 1) * DesPitch];
							auto Curr0L_S = DesS[(X + cPredDist) + (Y + 0) * DesPitch];
							auto Next1L_S = DesS[(X + cPredDist) + (Y + 1) * DesPitch];
							Temp = DeltaQubicWaveletPost1(Last2L_S, Last1L_S, Curr0L_S, Next1L_S);
						} else
						if (Y + 1 == SizeY0) {
							auto Last2L_S = DesS[(X + cPredDist) + (Y - 2) * DesPitch];
							auto Last1L_S = DesS[(X + cPredDist) + (Y - 1) * DesPitch];
							auto Curr0L_S = DesS[(X + cPredDist) + (Y + 0) * DesPitch];
							Temp = DeltaQubicWaveletPost2(Last2L_S, Last1L_S, Curr0L_S);
						}
						auto V = DesV[(X + cPredDist) + Y * DesPitch] - Temp;
						DesV[(X + cPredDist) + Y * DesPitch] = (tInt16)V;
					}
				}
				
				if (Y < 0 || X < 0 || (SizeX & SizeY) < 8) {
					continue;
				}
				
				auto TempH = (tInt16)0;
				auto TempD = (tInt16)0;
				if (X == 0) {
					TempH = DeltaQubicWaveletPre2(Curr0_S, Next1_S, Next2_S);
					TempD = DeltaQubicWaveletPre2(Curr0_V, Next1_V, Next2_V);
				} else
				if (X == 1) {
					TempH = DeltaQubicWaveletPre1(Last1_S, Curr0_S, Next1_S, Next2_S);
					TempD = DeltaQubicWaveletPre1(Last1_V, Curr0_V, Next1_V, Next2_V);
				} else
				if (X + 2 < SizeX1) {
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
					auto D = DesD[X + Y * DesPitch] - ((TempD + cDFactor) >> cDFactor);
					DesD[X + Y * DesPitch] = (tInt16)D;
				}
			}
		}
	}

	//=============================================================================
	static
	void ComposeWithQubPred(
		tLevel* SrcLevel,
		tLevel* DesLevel
	//=============================================================================
	) {
		MEASURE_BLOCK(ComposeWithQubPred);
		
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
		
		auto Last2_S = (tInt16)0;
		auto Last2_V = (tInt16)0;
		auto Last1_S = (tInt16)0;
		auto Last1_V = (tInt16)0;
		auto Curr0_S = (tInt16)0;
		auto Curr0_V = (tInt16)0;
		auto Next1_S = (tInt16)0;
		auto Next1_V = (tInt16)0;
		auto Next2_S = (tInt16)0;
		auto Next2_V = (tInt16)0;
		
		for (auto Y = 0; Y < SizeY0; Y++) {
			auto Y2 = Y << 1;
			for (auto X = -cPredDist; X < SizeX0; X++) {
				Last2_S = Last1_S;
				Last2_V = Last1_V;
				Last1_S = Curr0_S;
				Last1_V = Curr0_V;
				Curr0_S = Next1_S;
				Curr0_V = Next1_V;
				Next1_S = Next2_S;
				Next1_V = Next2_V;
				if (X + cPredDist < SizeX0 && (SizeX & SizeY) >= 8) {
					Next2_S = SrcS[(X + cPredDist) + Y * SrcPitch];
					
					if (Y < SizeY1) {
						auto Temp = (tInt16)0;
						if (Y == 0) {
							auto Curr0L_S = SrcS[(X + cPredDist) + (Y + 0) * SrcPitch];
							auto Next1L_S = SrcS[(X + cPredDist) + (Y + 1) * SrcPitch];
							auto Next2L_S = SrcS[(X + cPredDist) + (Y + 2) * SrcPitch];
							Temp = DeltaQubicWaveletPre2(Curr0L_S, Next1L_S, Next2L_S);
						} else
						if (Y == 1) {
							auto Last1L_S = SrcS[(X + cPredDist) + (Y - 1) * SrcPitch];
							auto Curr0L_S = SrcS[(X + cPredDist) + (Y + 0) * SrcPitch];
							auto Next1L_S = SrcS[(X + cPredDist) + (Y + 1) * SrcPitch];
							auto Next2L_S = SrcS[(X + cPredDist) + (Y + 2) * SrcPitch];
							Temp = DeltaQubicWaveletPre1(Last1L_S, Curr0L_S, Next1L_S, Next2L_S);
						} else
						if (Y + 2 < SizeY0) {
							auto Last2L_S = SrcS[(X + cPredDist) + (Y - 2) * SrcPitch];
							auto Last1L_S = SrcS[(X + cPredDist) + (Y - 1) * SrcPitch];
							auto Curr0L_S = SrcS[(X + cPredDist) + (Y + 0) * SrcPitch];
							auto Next1L_S = SrcS[(X + cPredDist) + (Y + 1) * SrcPitch];
							auto Next2L_S = SrcS[(X + cPredDist) + (Y + 2) * SrcPitch];
							Temp = DeltaQubicWavelet(Last2L_S, Last1L_S, Curr0L_S, Next1L_S, Next2L_S);
						} else
						if (Y + 2 == SizeY0) {
							auto Last2L_S = SrcS[(X + cPredDist) + (Y - 2) * SrcPitch];
							auto Last1L_S = SrcS[(X + cPredDist) + (Y - 1) * SrcPitch];
							auto Curr0L_S = SrcS[(X + cPredDist) + (Y + 0) * SrcPitch];
							auto Next1L_S = SrcS[(X + cPredDist) + (Y + 1) * SrcPitch];
							Temp = DeltaQubicWaveletPost1(Last2L_S, Last1L_S, Curr0L_S, Next1L_S);
						} else
						if (Y + 1 == SizeY0) {
							auto Last2L_S = SrcS[(X + cPredDist) + (Y - 2) * SrcPitch];
							auto Last1L_S = SrcS[(X + cPredDist) + (Y - 1) * SrcPitch];
							auto Curr0L_S = SrcS[(X + cPredDist) + (Y + 0) * SrcPitch];
							Temp = DeltaQubicWaveletPost2(Last2L_S, Last1L_S, Curr0L_S);
						}
						SrcV[(X + cPredDist) + Y * SrcPitch] += (tInt16)Temp;
					}
					Next2_V = SrcV[(X + cPredDist) + (Y + 0) * SrcPitch];
				}
				
				if (X >= 0 && (SizeX & SizeY) >= 8) {
					auto TempH = (tInt16)0;
					auto TempD = (tInt16)0;
					if (X == 0) {
						TempH = DeltaQubicWaveletPre2(Curr0_S, Next1_S, Next2_S);
						TempD = DeltaQubicWaveletPre2(Curr0_V, Next1_V, Next2_V);
					} else
					if (X == 1) {
						TempH = DeltaQubicWaveletPre1(Last1_S, Curr0_S, Next1_S, Next2_S);
						TempD = DeltaQubicWaveletPre1(Last1_V, Curr0_V, Next1_V, Next2_V);
					} else
					if (X + 2< SizeX1) {
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
						SrcD[X + Y * SrcPitch] += (tInt16)(TempD + cDFactor) >> cDFactor;
					}
				}
				
				if (X >= 0) {
					auto X2 = (tNat16)(X << 1);
					
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
							ASSERT(
								Res.A >= 0 &&
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
	}
	
#	ifdef TEST
		//=============================================================================
		static inline
		void TestSplitAndCompose(
			tArray<tNat8> Mem
		//=============================================================================
		) {
			const auto cSizeMax = (tNat16)128;
			const auto cDeltaPitch = (tNat16)16;
			
			auto Big = MemAlloc<tInt16>(Mem, (cSizeMax + cDeltaPitch) * cSizeMax);
			auto S = MemAlloc<tInt16>(Mem, (cSizeMax/2 + cDeltaPitch) * cSizeMax/2);
			auto H = MemAlloc<tInt16>(Mem, (cSizeMax/2 + cDeltaPitch) * cSizeMax/2);
			auto V = MemAlloc<tInt16>(Mem, (cSizeMax/2 + cDeltaPitch) * cSizeMax/2);
			auto D = MemAlloc<tInt16>(Mem, (cSizeMax/2 + cDeltaPitch) * cSizeMax/2);
			
			tLevel SmallLevel;
			SmallLevel.SizeXLeft   = cSizeMax/2;
			SmallLevel.SizeXRight  = cSizeMax/2;
			SmallLevel.SizeYTop    = cSizeMax/2;
			SmallLevel.SizeYBottom = cSizeMax/2;
			SmallLevel.S = S;
			SmallLevel.H = H;
			SmallLevel.V = V;
			SmallLevel.D = D;
			SmallLevel.Pitch = cSizeMax/2 + cDeltaPitch;
			
			tLevel BigLevel;
			BigLevel.SizeXLeft   = cSizeMax;
			BigLevel.SizeXRight  = cSizeMax;
			BigLevel.SizeYTop    = cSizeMax;
			BigLevel.SizeYBottom = cSizeMax;
			BigLevel.S = Big;
			BigLevel.Pitch = cSizeMax + cDeltaPitch;
			
			auto CalcPixel = [](tNat32 X, tNat32 Y, tNat32 I) { return (tInt16)(((X+257) * (Y+263) * I) % 251); };
			
			for (auto I = (tNat32)500; I -->0; ) {
				for (auto Y = (tNat32)0; Y < cSizeMax; Y += 1) {
					for (auto X = (tNat32)0; X < cSizeMax; X += 1) {
						Big[X + Y*cSizeMax] = CalcPixel(X, Y, I);
					}
				}
				
				for (auto Size = (tNat32)2; Size <= cSizeMax; Size <<= 1) {
					SplitWithQubPred  (&BigLevel, &SmallLevel);
					ComposeWithQubPred(&SmallLevel, &BigLevel);
					
					for (auto Y = (tNat16)0; Y < cSizeMax; Y += 1) {
						for (auto X = (tNat16)0; X < cSizeMax; X += 1) {
							auto A = Big[X + Y*cSizeMax];
							auto B = (tInt16)CalcPixel(X, Y, I);
							ASSERT(A == B);
						}
					}
				}
			}
		}
		
		//=============================================================================
		static inline
		void Test(
			tArray<tNat8> Mem
		//=============================================================================
		) {
			TestSplitAndCompose(Mem);
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
		tNat32 SizeX;
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
	tNat8 GetShift(
		tNat32 Mask
	//=============================================================================
	) {
		for (auto Shift = (tNat8)0; Shift < 32; Shift += 1) {
			if ((Mask & 1) != 0) {
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
//		MEASURE_BLOCK(GetLayerFromBmp16)
		
		auto MaskR = (tNat32)0x7C00;
		auto MaskG = (tNat32)0x03E0;
		auto MaskB = (tNat32)0x001F;
		if (BmpInfoHeader.Compression == 3) {
			GetColorMask(StreamIn, &MaskR, &MaskG, &MaskB);
		}
		auto MaskA = ~(MaskR | MaskG | MaskB) & 0x0000FFFF;
		
		auto ShiftA = GetShift(MaskA);
		auto ShiftR = GetShift(MaskR);
		auto ShiftG = GetShift(MaskG);
		auto ShiftB = GetShift(MaskB);
		
		auto SizeX = BmpInfoHeader.SizeX;
		auto RowSize = (tSize)(2*SizeX + 3) & (~3);
		auto Row = (tNat16*)malloc(RowSize);
		for (auto Y = (tNat16)0; fread(Row, RowSize, 1, StreamIn); Y += 1) {
			for (auto X = (tNat16)0; X < SizeX; X += 1) {
				auto J = Y*SizeX + X;
				LayerA->S[J] = (Row[X] & MaskA) >> ShiftA;
				
				auto R = (tInt16)((Row[X] & MaskR) >> ShiftR);
				auto G = (tInt16)((Row[X] & MaskG) >> ShiftG);
				auto B = (tInt16)((Row[X] & MaskB) >> ShiftB);
				
				auto RB = (tInt16)(R + B);
				
				LayerY->S[J]  = ((G<<1) + RB + 2) >> 2;
				LayerCg->S[J] = ((G<<1) - RB + 1) >> 1;
				LayerCo->S[J] = B - R;
			}
		}
		free(Row);
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
		MEASURE_BLOCK(GetLayerFromBmp24);
		
		auto SizeX = (tNat16)BmpInfoHeader.SizeX;
		auto Pitch = (tNat16)LayerA->Pitch;
		auto RowSize = (tSize)((3*SizeX+ 3) & (~3));
		auto Row = (tNat8*)malloc(RowSize);
		for (auto Y = (tNat32)0; fread(Row, RowSize, 1, StreamIn); Y += 1) {
			auto I = (tNat32)0;
			for (auto X = (tNat32)0; X < SizeX; X += 1) {
				auto J = Y*Pitch + X;
				LayerA->S[J] = 255;
				
				auto B = (tInt16)Row[I];
				I += 1;
				auto G = (tInt16)Row[I];
				I += 1;
				auto R = (tInt16)Row[I];
				I += 1;
				
				auto RB = (tInt16)(R + B);
				
				LayerY->S[J]  = ((G<<1) + RB + 2) >> 2;
				LayerCg->S[J] = ((G<<1) - RB + 1) >> 1;
				LayerCo->S[J] = B - R;
			}
		}
		free(Row);
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
		MEASURE_BLOCK(GetLayerFromBmp32);
		
		auto MaskR = (tNat32)0x00FF0000;
		auto MaskG = (tNat32)0x0000FF00;
		auto MaskB = (tNat32)0x000000FF;
		if (BmpInfoHeader.Compression == 3) {
			GetColorMask(StreamIn, &MaskR, &MaskG, &MaskB);
		}
		auto MaskA = ~(MaskR | MaskG | MaskB);
		
		auto ShiftA = GetShift(MaskA);
		auto ShiftR = GetShift(MaskR);
		auto ShiftG = GetShift(MaskG);
		auto ShiftB = GetShift(MaskB);
		
		auto SizeX = BmpInfoHeader.SizeX;
		auto RowSize = 4 * SizeX;
		auto Row = (tNat32*)malloc(RowSize);
		for (auto Y = (tNat32)0; fread(Row, RowSize, 1, StreamIn); Y += 1) {
			for (auto X = (tNat32)0; X < SizeX; X += 1) {
				auto J = (tNat32)(Y*SizeX + X);
				LayerA->S[J] = (tInt16)(Row[X] & MaskA) >> ShiftA;
				
				auto R = (tInt16)((Row[X] & MaskR) >> ShiftR);
				auto G = (tInt16)((Row[X] & MaskG) >> ShiftG);
				auto B = (tInt16)((Row[X] & MaskB) >> ShiftB);
				
				auto RB = (tInt16)(R + B);
				
				LayerY->S[J]  = ((G<<1) + RB + 2) >> 2;
				LayerCg->S[J] = ((G<<1) - RB + 1) >> 1;
				LayerCo->S[J] = B - R;
			}
		}
		free(Row);
	}
	
	//=============================================================================
	static
	void PutBmp32Header(
		tStream* StreamOut,
		tNat32   SizeX,
		tNat32   SizeY
	//=============================================================================
	) {
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
		
		tNat32 Mask[3] = {
			0x00FF0000u, // R
			0x0000FF00u, // G
			0x000000FFu  // B
		};
		if (!fwrite(&Mask, sizeof(Mask), 1, StreamOut)) {
			fprintf(gLogStream, "ERROR: fail writing to output!!!");
			exit(-1);
		}
	}
	
	//=============================================================================
	static
	void PutLayerToBmp32(
		tStream*       StreamOut,
		tNat32         SizeX,
		tNat32         SizeY,
		Layer::tLevel* LayerA,
		Layer::tLevel* LayerY,
		Layer::tLevel* LayerCg,
		Layer::tLevel* LayerCo
	//=============================================================================
	) {
		MEASURE_BLOCK(PutLayerToBmp32);
		
		auto MaskR = 0x00FF0000u;
		auto MaskG = 0x0000FF00u;
		auto MaskB = 0x000000FFu;
		auto MaskA = ~(MaskR | MaskG | MaskB);
		
		auto ShiftA = GetShift(MaskA);
		auto ShiftR = GetShift(MaskR);
		auto ShiftG = GetShift(MaskG);
		auto ShiftB = GetShift(MaskB);
		
		auto Pitch = LayerA->Pitch;
		auto RowSize = (tSize)(4*SizeX);
		auto Row = (tNat32*)malloc(RowSize);
		ASSERT(SizeY > 0); // TODO: allow negative SizeY
		for (auto Y = (tNat32)0; Y < (tNat32)SizeY; Y += 1) {
			for (auto X = (tNat32)0; X < SizeX; X += 1) {
				auto J = (tNat32)(Y*Pitch + X);
				
				auto A  = (tInt16)LayerA->S[J];
				auto Y_ = (tInt16)LayerY->S[J];
				auto Cg = (tInt16)LayerCg->S[J];
				auto Co = (tInt16)LayerCo->S[J];
				
				auto RB = (tInt16)((Y_<<1) - Cg);
				
				auto G = (tInt16)(((Y_<<1) + Cg) >> 1);
				auto R = (tInt16)((RB - Co + 1) >> 1);
				auto B = (tInt16)((RB + Co + 1) >> 1);
				
				auto Value = (
					Min<tInt16>(Max<tInt16>(0, R), 255) << ShiftR |
					Min<tInt16>(Max<tInt16>(0, G), 255) << ShiftG |
					Min<tInt16>(Max<tInt16>(0, B), 255) << ShiftB |
					Min<tInt16>(Max<tInt16>(0, A), 255) << ShiftA
				);
				Row[X] = Value;
			}
			fwrite(Row, RowSize, 1, StreamOut);
		}
		free(Row);
	}
}

using iCurve = struct {
	void* Env;
	tFunc3<void*, tNat32, tNat32, void>* Init;
	tFunc1<void*, void>*   Next;
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
		MemZero(State, 1);
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
		MemZero(State, 1);
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
		tStack Stack = 0;
		tNat32 X     = 0;
		tNat32 Y     = 0;
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
		End   = 0,
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
		static tStepType Matrix[] = { D, A, A, B };
		
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
		for (auto Temp = SizeX | SizeY; Temp != 0; Temp >>= 1) {
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
//		MEASURE_BLOCK(HilbertCurve_Next);
		
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
				ASSERT(cFalse);
			}
		}
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
			tArray<tNat8> Mem
		//=============================================================================
		) {
			const auto cBits = (tNat8)10;
			const auto cCount = (tNat32)(1 << (2*cBits));
			
			auto Array = MemAlloc<tBool>(Mem, cCount);
			MemZero(Array.Values, cCount);
			
			tState HilbertState;
			MemZero(&HilbertState, 1);
			HilbertState.Stack = HilbertCurve::New(cBits);
			for (auto Steps = cCount; Steps --> 0; ) {
				auto Pos = HilbertState.X | (HilbertState.Y << cBits);
				ASSERT(!Array[Pos]);
				Array[Pos] = cTrue;
				HilbertCurve::Next(&HilbertState);
			}
			
			for (auto I = cCount; I --> 0; ) {
				ASSERT(Array[I]);
			}
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
	ByteStream::tWriter* Stream,
	tArray<tInt16>       DataPtr,
	tNat32               SizeX,
	tNat32               SizeY,
	tNat32               Pitch,
	tNat8                Compression
//=============================================================================
) {
	MEASURE_BLOCK(WriteLayer);
	
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
	
	auto CompressionOffset = (tInt16)((1 << Compression) >> 1);
	
	tNat32 HistValues[1<<15];
	MemZero(HistValues, 1<<15);
	tNat32 HistZeros[1<<15];
	MemZero(HistZeros, 1<<15);
	
	for (auto I = (tNat32)0; I < SizeX*SizeY; ) {
		auto X = GetX(&Curve);
		auto Y = GetY(&Curve);
		if (X < SizeX && Y < SizeY) {
			I += 1;
			auto Value = (DataPtr[X + Y*Pitch] + CompressionOffset) >> Compression;
			if (Value == 0) {
				Zeros += 1;
			} else {
				if (Zeros > 0) {
					HistValues[0] += Zeros;
					HistZeros[UsedBits(Zeros)] += 1;
				}
				HistValues[Abs(Value)] += 1;
				
				if (Zeros > 0) {
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
	
	if (Zeros > 0) {
		HistValues[0] += Zeros;
		HistZeros[Min<tNat32>(Zeros, (1<<15) - 1)] += 1;
		if (Zeros == 1) {
			WriteNat(&BitWriterInterface, 0);
		} else {
			WriteNat(&BitWriterInterface, 1);
			WriteNat(&BitWriterInterface, Zeros - 2);
		}
	}
	
	Close(&BitWriterInterface);
	
#	if 0
	{
		MEASURE_BLOCK(WriteLayer_ShowHistogram);
		
		fprintf(gLogStream, "\n\n");
		fprintf(gLogStream, "%d x %d\n", SizeX, SizeY);
		for (auto Value = (tNat16)0; Value < 1<<15; Value += 1) {
			auto Count = HistValues[Value];
			if (Count > 0) {
				fprintf(gLogStream, "V %d %d\n", Value, Count);
			}
		}
		fprintf(gLogStream, "--\n");
		for (auto Repeates = (tNat32)0; Repeates < 1<<15; Repeates += 1) {
			auto Count = HistZeros[Repeates];
			if (Count > 0) {
				fprintf(gLogStream, "Z %d %d\n", Repeates, Count);
			}
		}
	}
#	endif
}

//=============================================================================
static
void ReadLayer(
	ByteStream::tReader* Stream,
	tArray<tInt16>       DataPtr,
	tNat32               SizeX,
	tNat32               SizeY,
	tNat32               Pitch
//=============================================================================
) {
	MEASURE_BLOCK(ReadLayer);
	
	BitStream::tReader BitStream;
	BitStream::Init(&BitStream, Stream);
	
	auto BitReaderInterface = GetInterface(&BitStream);
	
//	ArithmeticBitStream::tReader ABitReader;
//	ArithmeticBitStream::Init(&ABitReader, &BitReaderInterface_);
//	
//	auto BitReaderInterface = GetInterface(&ABitReader);
	
	auto Compression = (tNat8)ReadNat(&BitReaderInterface);
	
	auto Levels = (tNat32)0;
	for (auto Size = Max(SizeX, SizeY); Size > 0; Size >>= 1) {
		Levels += 1;
	}
	auto Zeros = (tNat64)0;
	auto WasLastZero = cFalse;
	tCurveState CurveState;
	auto Curve = GetInterface(&CurveState);
	Init(&Curve, SizeX, SizeY);
	
	for (auto I = (tNat32)0; I < SizeX * SizeY; ) {
		auto X = GetX(&Curve);
		auto Y = GetY(&Curve);
		if (X < SizeX && Y < SizeY) {
			I += 1;
			tInt8 Sign;
			tNat64 Value;
			if (Zeros > 0) {
				Zeros -= 1;
				Sign = 1;
				Value = 0;
			} else if (WasLastZero) {
				Value = ReadNat(&BitReaderInterface);
				Sign = ((Value & 1) == 0) ? 1 : -1;
				Value >>= 1;
				Value += 1;
			} else {
				Value = ReadNat(&BitReaderInterface);
				Sign = ((Value & 1) == 0) ? 1 : -1;
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
}

//=============================================================================
void EnCode(
	tStream*             StreamIn,
	ByteStream::tWriter* StreamOut,
	tNat32               Compression,
	tArray<tNat8>        Mem
//=============================================================================
) {
	MEASURE_BLOCK(EnCode);
	
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
	
	Layer::tLevel Layers[4][32];
	auto LayerCount = (tNat16)4;
	
	auto InitLevel = (tNat32)0;
	{
		for (auto Layer = LayerCount; Layer --> 0; ) {
			InitLevel = Layer::Init(Layers[Layer], (tNat16)BmpInfoHeader.SizeX, (tNat16)BmpInfoHeader.SizeY, Mem);
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
	for (auto Layer = LayerCount; Layer --> 0; ) {
		for (auto Level = InitLevel; Level --> 0; ) {
			Layer::SplitWithQubPred(&Layers[Layer][Level+1], &Layers[Layer][Level]);
		}
	}
	
#	ifdef DEBUG
	{
		for (auto Layer = LayerCount; Layer --> 0; ) {
			for (auto Level = InitLevel - 1; Level --> 0; ) {
				for (auto Part = 0; Part < 3; Part += 1) {
					auto LayerPtr = &Layers[Layer][Level+1];
					auto Pitch = LayerPtr->Pitch;
					tNat16 SizeX = 0;
					tNat16 SizeY = 0;
					tArray<tInt16> Data = {};
					auto Mode = '_';
					switch (Part) {
						case 0: {
							Mode = 'H';
							Data = LayerPtr->H;
							SizeX = LayerPtr->SizeXRight;
							SizeY = LayerPtr->SizeYTop;
						} break;
						
						case 1: {
							Mode = 'V';
							Data = LayerPtr->V;
							SizeX = LayerPtr->SizeXLeft;
							SizeY = LayerPtr->SizeYBottom;
						} break;
						
						case 2: {
							Mode = 'D';
							Data = LayerPtr->D;
							SizeX = LayerPtr->SizeXRight;
							SizeY = LayerPtr->SizeYBottom;
						} break;
						
						default: ASSERT(false);
					}
					char str[] = {
						(char)('0' + (Level / 10) % 10),
						(char)('0' + (Level /  1) % 10),
						'_',
						(char)('0' + (Layer / 10) % 10),
						(char)('0' + (Layer /  1) % 10),
						'_',
						Mode,
						'.',
						'b',
						'm',
						'p',
						'\0'
					};
					
					auto StreamOut = fopen(str, "wb");
					
					BmpHelper::PutBmp32Header(StreamOut, SizeX, SizeY);
					
					auto RowSize = 4 * SizeY;
					auto Row = (tNat32*)malloc(RowSize);
					for (auto Y = (tNat32)0; Y < (tNat32)SizeY; Y += 1) {
						for (auto X = (tNat32)0; X < SizeX; X += 1) {
							auto Value = Data[X + Y*Pitch] << 0;
							
							auto R = (tInt16)(-Value);
							auto G = (tInt16)(Value);
							auto B = (tInt16)((Abs(Value) > 0xFF) ? 0xFF : 0);
							
							Row[X] = (
								((tInt32)((Min<tInt16>(Max<tInt16>(0, R), 0xFF)) << 16)) |
								((tInt32)((Min<tInt16>(Max<tInt16>(0, G), 0xFF)) <<  8)) |
								((tInt32)((Min<tInt16>(Max<tInt16>(0, B), 0xFF)) <<  0))
							);
						}
						fwrite(Row, RowSize, 1, StreamOut);
					}
					
					auto fclose(StreamOut);
				}
			}
		}
	}
#	endif
	
	{ // Output
		auto InitLevelPtr = &Layers[0][InitLevel];
		auto SizeX = InitLevelPtr->SizeXLeft;
		auto SizeY = InitLevelPtr->SizeYTop;

		BitStream::tWriter BitWriter;
		BitStream::Init(&BitWriter, StreamOut);
		auto BitWriterInterface = GetInterface(&BitWriter);
		
		WriteNat(&BitWriterInterface, SizeX);
		WriteNat(&BitWriterInterface, SizeY);
		
		for (auto Layer = LayerCount; Layer --> 0; ) {
			auto Value = Layers[Layer][0].S[0];
			WriteBit(&BitWriterInterface, (Value < 0) ? 1 : 0);
			WriteNat(&BitWriterInterface, Abs(Value));
			
			ASSERT(Layers[Layer][0].SizeXLeft == 1);
			ASSERT(Layers[Layer][0].SizeYTop == 1);
		}
		
		Close(&BitWriterInterface);
		
		for (auto Level = (tNat32)0; Level < InitLevel; Level += 1) {
			auto SizeXLeft   = Layers[0][Level].SizeXLeft;
			auto SizeXRight  = Layers[0][Level].SizeXRight;
			auto SizeYTop    = Layers[0][Level].SizeYTop;
			auto SizeYBottom = Layers[0][Level].SizeYBottom;
			auto Pitch       = Layers[0][Level].Pitch;
			
			ASSERT(SizeXLeft >= SizeXRight);
			ASSERT(SizeYTop >= SizeYBottom);
			
			ASSERT(SizeXLeft + SizeXRight == Layers[0][Level+1].SizeXLeft);
			ASSERT(SizeYTop + SizeYBottom == Layers[0][Level+1].SizeYTop);
			
			for (auto Layer = (tNat32)0; Layer < LayerCount; Layer += 1) {
				auto LayerTemp = Layers[Layer][Level];
				
				auto Compression_ = (tNat8)Max<tInt32>(0, Level - InitLevel + Compression - 2*(Layer == 1));
				
				WriteLayer(StreamOut, LayerTemp.H, SizeXRight, SizeYTop,    Pitch, Compression_);
				WriteLayer(StreamOut, LayerTemp.V, SizeXLeft,  SizeYBottom, Pitch, Compression_);
				WriteLayer(StreamOut, LayerTemp.D, SizeXRight, SizeYBottom, Pitch, Compression_);
			}
		}
	}
}

//=============================================================================
void DeCode(
	ByteStream::tReader* StreamIn,
	tStream*             StreamOut,
	tArray<tNat8>        Mem
//=============================================================================
) {
	MEASURE_BLOCK(DeCode);
	
	BitStream::tReader BitReader;
	BitStream::Init(&BitReader, StreamIn);
	auto BitReaderInterface = GetInterface(&BitReader);
	
	auto SizeX = (tNat16)ReadNat(&BitReaderInterface);
	auto SizeY = (tNat16)ReadNat(&BitReaderInterface);
	
	Layer::tLevel Layers[4][32];
	auto LayerCount = (tNat8)4;
	auto InitLevel = (tNat8)0;
	{
		auto BufferSize = 4 * SizeX * Max<tNat32>(SizeY, 8);
		auto Buffer = MemAlloc<tInt16>(Mem, LayerCount * BufferSize * sizeof(tInt16));
		
		{
			for (auto Layer = LayerCount; Layer --> 0; ) {
				InitLevel = Layer::Init(&Layers[Layer][0], SizeX, SizeY, Mem);
			}
		}
		
		for (auto Layer = LayerCount; Layer --> 0; ) {
			auto Sign = (tInt16)(ReadBit(&BitReaderInterface) == 1 ? -1 : 1);
			Layers[Layer][0].S[0] = (tInt16)(Sign * ReadNat(&BitReaderInterface));
			
			ASSERT(Layers[Layer][0].SizeXLeft == 1);
			ASSERT(Layers[Layer][0].SizeYTop == 1);
		}
		
		Close(&BitReaderInterface);
		
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
				
				ReadLayer(StreamIn, LayerTemp.H, SizeXRight, SizeYTop,    Pitch);
				ReadLayer(StreamIn, LayerTemp.V, SizeXLeft,  SizeYBottom, Pitch);
				ReadLayer(StreamIn, LayerTemp.D, SizeXRight, SizeYBottom, Pitch);
			}
		}
	}
	
	{ // Berechnen
		for (auto Layer = LayerCount; Layer --> 0; ) {
			for (auto Level = (tNat16)0; Level < InitLevel; Level += 1) {
				Layer::ComposeWithQubPred(&Layers[Layer][Level], &Layers[Layer][Level+1]);
			}
		}
	}
	
	{ // Output
		BmpHelper::PutBmp32Header(StreamOut, SizeX, SizeY);
		BmpHelper::PutLayerToBmp32(
			StreamOut,
			SizeX,
			SizeY,
			&Layers[0][InitLevel],
			&Layers[1][InitLevel],
			&Layers[2][InitLevel],
			&Layers[3][InitLevel]
		);
	}
}

#ifdef TEST
	//=============================================================================
	static inline
	void Test(
		tArray<tNat8> Mem
	//=============================================================================
	) {
		TestArray();
		HaarWavelet::Test();
		HilbertCurve::Test(Mem);
		FibonacciCode::Test();
		Huffman::Test();
		Layer::Test(Mem); 
	}
#endif

tNat8 gMem[1<<25];

//=============================================================================
int main(
	int      ArgCount,
	tChar8** Args
//=============================================================================
) {
	START_CLOCK(main);
	
	gLogStream = stderr;
	
	FibonacciCode::Init(FibonacciCode::gFibArray, _countof(FibonacciCode::gFibArray));
	
	tArray<tNat8> Mem = AsArray(gMem);
	
#	ifdef TEST
		if (ArgCount > 1 && strncmp(Args[1], "-t", 2) == 0) {
			
			Test(Mem);
			if (ArgCount == 2) {
				exit(0);
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
		if (!StreamOut) {
			fprintf(gLogStream, "ERROR: Can't open/create File '%s'!!!", FileInName);
			exit(-1);
		}
	}
	
	if (strncmp(Flag, "-e", 2) == 0) {
		BufferdStream::tWriter BufferdStreamOut;
		BufferdStream::Init(&BufferdStreamOut, StreamOut);
		auto x = ByteStream::Writer(
			&BufferdStreamOut,
			[] (void* Env, tNat8 Byte) { BufferdStream::WriteByte((BufferdStream::tWriter*)Env, Byte); },
			[] (void* Env) { BufferdStream::Close((BufferdStream::tWriter*)Env); }
		);
		EnCode(StreamIn, &x, atoi(&Flag[2]), Mem);
		BufferdStream::Close(&BufferdStreamOut);
	} else
	if (strncmp(Flag, "-d", 2) == 0) {
		BufferdStream::tReader BufferdStreamIn;
		BufferdStream::Init(&BufferdStreamIn, StreamIn);
		auto x = ByteStream::Reader(
			&BufferdStreamIn,
			[] (void* Env) { return BufferdStream::ReadByte((BufferdStream::tReader*)Env); },
			[] (void* Env) { BufferdStream::Close((BufferdStream::tReader*)Env); }
		);
		DeCode(&x, StreamOut, Mem);
		BufferdStream::Close(&BufferdStreamIn);
	} else {
		fprintf(gLogStream, "ERROR: unknown arg %s!!!", Flag);
		exit(-1);
	}
	
	STOP_CLOCK(main);
#	ifdef VS
		PRINT_CLOCKS(fopen("C:\\Projekte\\ImgCodec2015\\bin\\timer.txt", "w"));
#	else
		PRINT_CLOCKS(gLogStream);
#	endif
	return 0;
}
