package ckks_cipher

import (
	"fmt"
	"math"
	"strconv"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

type Ascon struct {
    *RtBCipher
	iv         	[]byte
    blockSize	int
	keySize	    int
	rounds		int
	// bitIndex   	[]*BitSet
    bitSbox	   	[]*BitSet
    sboxMonomialOrder	[]*BitSet
	allZeroIn	bool
}

func NewAscon(key_ []uint8, params_ ckks.Parameters, btpParams_ bootstrapping.Parameters, btpKey_ *bootstrapping.EvaluationKeys, encoder_ *ckks.Encoder, encryptor_ *rlwe.Encryptor, decryptor_ *rlwe.Decryptor, iv_ []byte) (*Ascon, error) {
    rtb, err := NewRtBCipher(key_, params_, btpParams_, btpKey_, encoder_, encryptor_, decryptor_)
    if err!= nil {
        return nil, err
    }
    aes := &Ascon{
        RtBCipher: 	rtb,
        blockSize:	128,
		keySize: 	128,
		rounds:		10,
        iv:         iv_,
		allZeroIn: 	true, //debug mode
    }

	return aes, err
}

func (aes *Ascon) DebugTest(ciphertexts []byte, bits int) ([]*rlwe.Ciphertext, error) {
	if aes.allZeroIn {
		bits = aes.params.MaxSlots() * aes.blockSize
	}
	numBlocks := int(math.Ceil(float64(bits) / float64(aes.blockSize)))

	iv := NewBitSet(aes.blockSize)
	iv.Set(0) // set iv all zero
	aes.EncryptKey()
	aes.EncryptInput(iv, numBlocks)
	state := aes.inputEncrypted

	fmt.Println("aes.totalLevel", aes.totalLevel)
	fmt.Println("aes.remain", aes.remainingLevel)
	for i := 0; i < len(state); i++ {
		if state[i] == nil { 
			strMsg := "state[" + strconv.Itoa(i) + "] is nil, exit!"	
			panic(strMsg) 
		}

		aes.DropLevel(state[i], aes.totalLevel-aes.remainingLevel)
	}
	aes.DebugPrint(state[0], "Level")
	// // aes.Evaluator.Add(state[0], state[0], state[0])
	// for i:=0;i<8;i++{
	// 	aes.DebugPrint(state[i], "before Bootstrap bits: \n")
	// }	
	// aes.aesRoundFunction(state[:], aes.keyEncrypted)
	Start := time.Now()

	ch := make(chan bool, 16)
	fmt.Println("len(state): ", len(state))
	for i:=0; i<16; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			aes.aesSubbyteLUT(evalCopy, state[i*8 : (i+1)*8])
			ch <- true
		}(i)
	}
	for i := 0; i < 16; i++ {
		<-ch
	}

	// aes.RoundFunction(evals, state, aes.keyEncrypted)
	elasp := time.Since(Start)
	fmt.Println("\n\n\n\nOne round time: ", elasp)
	// valuesTest := (*aes.encoder).DecodeComplex( (*aes.decryptor).DecryptNew(state[0]), aes.params.LogSlots())
	// fmt.Println("BootReEnc debug")
	// PrintVectorTrunc(valuesTest, 7, 3)
	for i:=0;i<8;i++{
		aes.DebugPrint(state[i], "0~after btp")
		
	}
	// for i:=0;i<8;i++{
	// 	aes.DebugPrint(state[i], "After One Round: \n")
	// }
	return state, nil
}

func (aes *Ascon) HEDecrypt(ciphertexts []uint8, bits int) []*rlwe.Ciphertext {
	if aes.allZeroIn {
		bits = aes.params.MaxSlots() * aes.blockSize
	}
	numBlock := int(math.Ceil(float64(bits) / float64(aes.blockSize)))
	iv := NewBitSet(aes.blockSize)
	iv.Set(0) // set iv all zero
	aes.EncryptKey()
	aes.EncryptInput(iv, numBlock)

	for i := 0; i < len(aes.inputEncrypted); i++ {
		if aes.inputEncrypted[i].Level() > aes.remainingLevel + 1  {
			aes.Evaluator.DropLevel(aes.inputEncrypted[i], aes.inputEncrypted[i].Level() - aes.remainingLevel -1 )
		}
	}

	startAES := time.Now()
	// AES encryption **********************************************
	state := aes.AddWhiteKey( aes.inputEncrypted, aes.keyEncrypted )
	// state := aes.inputEncrypted
	for i := 1; i < 10; i++ {
		fmt.Printf("round iterator : %d\n", i)
		aes.RoundFunction(state, aes.keyEncrypted)
	}
	fmt.Println("round iterator : last round")
	aes.LastRound( state, aes.keyEncrypted)
	// AES encryption **********************************************
	for i:=0;i<8;i++{
		str := "Sbox: " + strconv.Itoa(i)
		aes.DebugPrint(state[i], str)
	}

	endAES := time.Now()
	durationAES := endAES.Sub(startAES)
	fmt.Printf("代码执行时间：%d 秒 :: %d 毫秒\n", int(durationAES.Seconds()), int(durationAES.Milliseconds())%1000)
	// // Add cipher
	// // encode_ciphertext(ciphertexts, num_block);
	// for i := 0; i < len(state); i++ {
	// 	if (i%8 == 1) || (i%8 == 2) || (i%8 == 4) || (i%8 == 5) {
	// 		NOT( evals[0], state[i], state[i])
	// 	}
	// }

	// ch := make(chan int, len(state) )
	// for i := 0; i < len(state); i++ {
	// 	go func(i int) {
	// 		XOR(evals[i], state[i], aes.encodeCipher[i], state[i])
	// 		ch <- i
	// 	}(i)
	// }
	// for i := 0; i < len(state); i++ {
	// 	<-ch
	// }

	return state
}


func (aes *Ascon) HEDecryptParam15(ciphertexts []uint8, bits int) []*rlwe.Ciphertext {
	if aes.allZeroIn {
		bits = aes.params.MaxSlots() * aes.blockSize
	}
	numBlock := int(math.Ceil(float64(bits) / float64(aes.blockSize)))
	iv := NewBitSet(aes.blockSize)
	iv.Set(0) // set iv all zero
	aes.EncryptKey()
	aes.EncryptInput(iv, numBlock)

	for i := 0; i < len(aes.inputEncrypted); i++ {
		if aes.inputEncrypted[i].Level() > aes.remainingLevel + 1 {
			aes.Evaluator.DropLevel(aes.inputEncrypted[i], aes.inputEncrypted[i].Level() - aes.remainingLevel - 1 )
		}
	}

	startAES := time.Now()
	state := aes.AddWhiteKey( aes.inputEncrypted, aes.keyEncrypted )
	// state := aes.inputEncrypted
	// AES encryption **********************************************
	for i := 1; i < 5; i++ {
		fmt.Printf("round iterator : %d\n", i)
		aes.RoundFunctionParam15(state, aes.keyEncrypted)
	}
	fmt.Println("round iterator : last round")
	aes.LastRoundParam15( state, aes.keyEncrypted)
	// AES encryption **********************************************
	for i:=0;i<8;i++{
		str := "Sbox: " + strconv.Itoa(i)
		aes.DebugPrint(state[i], str)
	}

	endAES := time.Now()
	durationAES := endAES.Sub(startAES)
	fmt.Printf("代码执行时间：%d 秒 :: %d 毫秒\n", int(durationAES.Seconds()), int(durationAES.Milliseconds())%1000)
	// // Add cipher
	// // encode_ciphertext(ciphertexts, num_block);
	// for i := 0; i < len(state); i++ {
	// 	if (i%8 == 1) || (i%8 == 2) || (i%8 == 4) || (i%8 == 5) {
	// 		NOT( evals[0], state[i], state[i])
	// 	}
	// }

	// ch := make(chan int, len(state) )
	// for i := 0; i < len(state); i++ {
	// 	go func(i int) {
	// 		XOR(evals[i], state[i], aes.encodeCipher[i], state[i])
	// 		ch <- i
	// 	}(i)
	// }
	// for i := 0; i < len(state); i++ {
	// 	<-ch
	// }

	return state
}

func (aes *Ascon) EncryptKey() {
    if aes.encoder == nil || aes.encryptor == nil {
        panic("encoder or encryptor is not initialized")
    }
    aes.keyEncrypted = make([]*rlwe.Ciphertext, aes.keySize)

    for i := 0; i < aes.keySize; i++ {
        if i >= len(aes.symmetricKey)*8 {
            panic("input symmetric key size is not match!")
        }
        bit := (aes.symmetricKey[i/8] >> uint(i%8)) & 1
        keys := make([]float64, aes.params.MaxSlots())
        for j := 0; j < aes.params.MaxSlots(); j++ {
            keys[j] = float64(bit)
        }
        // fmt.Printf("Encoding bit %d which is %d\n", i, bit)
		pt := ckks.NewPlaintext( aes.params, aes.totalLevel)
        aes.encoder.Encode(keys, pt)
        if pt == nil {
            panic("skPlain is nil after encoding")
        }
        // fmt.Println("Before encryption")
        aes.keyEncrypted[i], _ = aes.encryptor.EncryptNew(pt)		
    }
	if aes.keyEncrypted[0] == nil {
		panic("skBitEncrypted is nil after encryption")
	}
    // fmt.Println("Encryption completed")
}

func (aes *Ascon) EncryptInput(iv *BitSet, numBlock int) {
	aes.inputEncrypted = make([]*rlwe.Ciphertext, aes.blockSize)
	inputData := make([]*BitSet, aes.params.MaxSlots() )
	for i, block := range inputData{
		block = NewBitSet(aes.blockSize)
		inputData[i] = block
	}

	for i := range inputData {
		if aes.allZeroIn {
			inputData[i].Set(0)
		} else {
			inputData[i].Set( int(ctr(iv, uint64(i+1) ).ToULong()) )
		}
	}
	for i := 0; i < aes.blockSize; i++ {
		stateBatched := make([]float64, aes.params.MaxSlots())
		for j := 0; j < aes.params.MaxSlots(); j++ {
			stateBatched[j] = float64(inputData[j].bits[i])
		}
		pt := ckks.NewPlaintext( aes.params, aes.totalLevel)
		aes.encoder.Encode( stateBatched, pt)
		aes.inputEncrypted[i], _ = aes.encryptor.EncryptNew(pt)
	}
	if aes.inputEncrypted[0] == nil {panic("input is not stored in aesStruct")}
}

func (aes *Ascon) EncodeCiphertext(ciphertexts []uint8, numBlock int) {
	aes.encodeCipher = make([]*rlwe.Ciphertext, aes.blockSize)
	if numBlock < aes.params.MaxSlots() {
		fmt.Println("data is not full pack, fill with 0...")
	}
	encryptedData := make([]*BitSet, numBlock)
	for i, bit := range encryptedData {
		bit = NewBitSet(aes.blockSize)
		encryptedData[i] = bit 
	}

	for i := range encryptedData {
		if i < numBlock {
			if aes.allZeroIn {
				encryptedData[i].Set(0)
			} else {
				for k := 0; k < aes.blockSize && i*aes.blockSize+k < numBlock*aes.blockSize; k++ {
					ind := i*aes.blockSize + k
					bit := (ciphertexts[ind/8] >> uint(ind%8)) & 1
					encryptedData[i].bits[k] = uint8(bit)
				}
			}
		} else {
			encryptedData[i].Set(0)
		}
	}
	for i := 0; i < aes.blockSize; i++ {
		var data []float64
		for j := 0; j < aes.params.MaxSlots(); j++ {
			data = append(data, float64( encryptedData[j].bits[i] ) )
		}
		pt := ckks.NewPlaintext( aes.params, aes.remainingLevel)
		aes.encoder.Encode(data, pt)
		aes.encodeCipher[i], _ = aes.encryptor.EncryptNew(pt)
	}
	if aes.encodeCipher[0] == nil {panic("encodeCiphertext is nil")}

}

func (aes *Ascon) RoundFunction( state []*rlwe.Ciphertext, roundKey []*rlwe.Ciphertext) {
	// SubByte
	fmt.Printf("Chain index before sbox: %d, scale: %f\n", state[77].Level(), state[77].LogScale() )
	ch := make(chan bool, 16)
	for i := 0; i < 16; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			aes.aesSubbyteLUT(evalCopy, state[i*8 : (i+1)*8])
			ch <- true
		}(i)
	}
	for i := 0; i < 16; i++ {
		<-ch
	}
	// ShiftRow
	aes.ShiftRow(state)
	// MixColumn
	aes.FreeMixColumn( state )
	fmt.Printf("MixColumn Chain: %d, scale: %f\n", state[0].Level(), state[0].LogScale() )
	// AddRoundKey
	aes.FreeAddRoundKey( state, roundKey )
	fmt.Printf("AddRoundKey Chain: %d, scale: %f\n", state[0].Level(), state[0].LogScale() )
	// Parallel processing for bootstrapping and cleaning tensor
	ch = make(chan bool, 64)
	for i:=0; i<64; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			state[i], state[64 + i], _ = evalCopy.BootstrapCmplxThenDivide(state[i], state[64 + i])
			if i == 0 {
				aes.DebugPrint(state[i], "BTS precise: ")
			}
			ch <- true
		}(i)
	}
	for i := 0; i < 64; i++ {
		<-ch
	}
}

func (aes *Ascon) LastRound( state []*rlwe.Ciphertext, roundKey []*rlwe.Ciphertext) {
	fmt.Printf("Chain index before sbox: %d, scale: %f\n", state[0].Level(), state[0].LogScale() )
	ch := make(chan bool, 16)
	for i := 0; i < 16; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			aes.aesSubbyteLUT(evalCopy, state[i*8 : (i+1)*8])
			ch <- true
		}(i)
	}
	for i := 0; i < 16; i++ {
		<-ch
	}
	// ShiftRow
	aes.ShiftRow(state)
	// AddRoundKey
	aes.FreeAddRoundKey( state, roundKey )

	ch = make(chan bool, 64)
	for i:=0; i<64; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			state[i], state[64 + i], _ = evalCopy.BootstrapCmplxThenDivide(state[i], state[64 + i])	
			if i == 0 {
				aes.DebugPrint(state[i], "BTS precise: ")
			}
			ch <- true
		}(i)
	}
	for i := 0; i < 64; i++ {
		<-ch
	}

}

func (aes *Ascon) RoundFunctionParam15( state []*rlwe.Ciphertext, roundKey []*rlwe.Ciphertext) {
	// SubByte first round
	ch := make(chan bool, 16)
	for i := 0; i < 16; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			aes.aesSubbyteLUT(evalCopy, state[i*8 : (i+1)*8])
			ch <- true
		}(i)
	}
	for i := 0; i < 16; i++ {
		<-ch
	}
	aes.DebugPrint(state[0], "After Sbox: ")
	// ShiftRow
	aes.ShiftRow(state)
	// MixColumn
	aes.MixColumn( state )
	fmt.Printf("First Round MixColumn Chain: %d, scale: %f\n", state[0].Level(), state[0].LogScale() )
	// AddRoundKey
	aes.AddRoundKey( state, roundKey )
	fmt.Printf("First Round AddRoundKey Chain: %d, scale: %f\n", state[0].Level(), state[0].LogScale() )
	// SubByte free XOR round
	for i := 0; i < 16; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			aes.aesSubbyteLUT(evalCopy, state[i*8 : (i+1)*8])
			ch <- true
		}(i)
	}
	for i := 0; i < 16; i++ {
		<-ch
	}
	aes.DebugPrint(state[0], "Second After Sbox: ")
	// ch = make(chan bool, 128)
	// for i:=0; i<128; i++ {
	// 	go func(i int) {
	// 		evalCopy := aes.Evaluator.ShallowCopy()
	// 		CleanReal(evalCopy, state[i])
	// 		ch <- true
	// 	}(i)
	// }
	// for i := 0; i < 128; i++ {
	// 	<-ch
	// }
	// aes.DebugPrint(state[0], "Clean Precise: ")
	// ShiftRow
	aes.ShiftRow(state)
	// MixColumn
	aes.FreeMixColumn( state )
	// aes.DebugPrint(state[0], "Third After MixColumn: ")
	fmt.Printf("Free MixColumn Chain: %d, scale: %f\n", state[0].Level(), state[0].LogScale() )
	// AddRoundKey
	aes.FreeAddRoundKey( state, roundKey )
	aes.DebugPrint(state[0], "Third After FreeAddKey: ")
	// fmt.Printf("Free AddRoundKey Chain: %d, scale: %f\n", state[0].Level(), state[0].LogScale() )
	ch = make(chan bool, 64)
	for i:=0; i<64; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			state[i], state[64 + i], _ = evalCopy.BootstrapCmplxThenDivide(state[i], state[64 + i])
			ch <- true
		}(i)
	}
	for i := 0; i < 64; i++ {
		<-ch
	}
	aes.DebugPrint(state[0], "BTS precise: ")
}

func (aes *Ascon) LastRoundParam15( state []*rlwe.Ciphertext, roundKey []*rlwe.Ciphertext) {
	// SubByte first round
	ch := make(chan bool, 16)
	for i := 0; i < 16; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			aes.aesSubbyteLUT(evalCopy, state[i*8 : (i+1)*8])
			ch <- true
		}(i)
	}
	for i := 0; i < 16; i++ {
		<-ch
	}
	aes.DebugPrint(state[0], "After Sbox: ")
	// ShiftRow
	aes.ShiftRow(state)
	// MixColumn
	aes.MixColumn( state )
	fmt.Printf("First Round MixColumn Chain: %d, scale: %f\n", state[0].Level(), state[0].LogScale() )
	// AddRoundKey
	aes.AddRoundKey( state, roundKey )
	fmt.Printf("First Round AddRoundKey Chain: %d, scale: %f\n", state[0].Level(), state[0].LogScale() )
	// SubByte free XOR round
	for i := 0; i < 16; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			aes.aesSubbyteLUT(evalCopy, state[i*8 : (i+1)*8])
			ch <- true
		}(i)
	}
	for i := 0; i < 16; i++ {
		<-ch
	}
	// ch = make(chan bool, 128)
	// for i:=0; i<128; i++ {
	// 	go func(i int) {
	// 		evalCopy := aes.Evaluator.ShallowCopy()
	// 		CleanReal(evalCopy, state[i])
	// 		ch <- true
	// 	}(i)
	// }
	// for i := 0; i < 128; i++ {
	// 	<-ch
	// }
	// aes.DebugPrint(state[0], "Clean Precise: ")
	// ShiftRow
	aes.ShiftRow(state)
	// AddRoundKey
	aes.FreeAddRoundKey( state, roundKey )
	ch = make(chan bool, 64)
	for i:=0; i<64; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			state[i], state[64 + i], _ = evalCopy.BootstrapCmplxThenDivide(state[i], state[64 + i])
			ch <- true
		}(i)
	}
	for i := 0; i < 64; i++ {
		<-ch
	}
}

func (aes *Ascon) coefficientMultMonomial(eval *bootstrapping.Evaluator, mon []*rlwe.Ciphertext, coeffArr []int, pos int) ( ctOut *rlwe.Ciphertext ) {
    if len(mon)!= len(aes.sboxMonomialOrder) {
		panic("monomial size must equal to sbox_monomial_order!")
    }
	ctOut = mon[0].CopyNew()
	i := 0
    for i < len(mon){
        ind := int(aes.sboxMonomialOrder[i].ToULong()) - 1
        coeff := coeffArr[ind]
		if coeff == 0 {
			i++
			continue
		} 
		eval.Mul(mon[i], coeff, ctOut)
		i++
		break
    }

	for i < len(mon){
        ind := int(aes.sboxMonomialOrder[i].ToULong()) - 1
        coeff := coeffArr[ind]
		if coeff == 0 {
			i++
			continue
		} 
		tmp := ctOut.CopyNew()
		eval.Mul(mon[i], coeff, tmp)
		eval.Add(tmp, ctOut, ctOut)
		i++
    }
	eval.Add( ctOut, int( aes.bitSbox[0].bits[pos] ), ctOut )
    return 
}

func (aes *Ascon) aesSubbyteLUT(eval *bootstrapping.Evaluator, SBoxIn []*rlwe.Ciphertext) {
    // construct 8-bit val of the sbox
    if len(SBoxIn)!= 8 {
		panic("The input length of the Sbox is wrong (8bit)!!")
    }
    sboxMonomials, _ := LayeredCombine(eval, SBoxIn)

    SBoxIn[0] = aes.coefficientMultMonomial(eval, sboxMonomials, Sbox0[:], 0)
    SBoxIn[1] = aes.coefficientMultMonomial(eval, sboxMonomials, Sbox1[:], 1)
    SBoxIn[2] = aes.coefficientMultMonomial(eval, sboxMonomials, Sbox2[:], 2)
    SBoxIn[3] = aes.coefficientMultMonomial(eval, sboxMonomials, Sbox3[:], 3)
    SBoxIn[4] = aes.coefficientMultMonomial(eval, sboxMonomials, Sbox4[:], 4)
    SBoxIn[5] = aes.coefficientMultMonomial(eval, sboxMonomials, Sbox5[:], 5)
    SBoxIn[6] = aes.coefficientMultMonomial(eval, sboxMonomials, Sbox6[:], 6)
    SBoxIn[7] = aes.coefficientMultMonomial(eval, sboxMonomials, Sbox7[:], 7)
}


func (aes *Ascon) MixColumn( x []*rlwe.Ciphertext ) {
	x0, x1, x2, x3 := []*rlwe.Ciphertext{}, []*rlwe.Ciphertext{}, []*rlwe.Ciphertext{}, []*rlwe.Ciphertext{}
	
	for i := 0; i < 128; i++ {
		mod := i % 32
		if mod < 8 {
			x0 = append(x0, x[i])
		} else if mod >= 8 && mod < 16 {
			x1 = append(x1, x[i])
		} else if mod >= 16 && mod < 24 {
			x2 = append(x2, x[i])
		} else {
			x3 = append(x3, x[i])
		}
	}

	y0, y1, y2, y3 := make([]*rlwe.Ciphertext, len(x0)), make([]*rlwe.Ciphertext, len(x1)), make([]*rlwe.Ciphertext, len(x2)), make([]*rlwe.Ciphertext, len(x3))
	z0, z1, z2, z3 := make([]*rlwe.Ciphertext, len(x0)), make([]*rlwe.Ciphertext, len(x1)), make([]*rlwe.Ciphertext, len(x2)), make([]*rlwe.Ciphertext, len(x3))
	
	
	ch := make(chan bool, 32)
	for i := 0; i < 32; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			y0[i] = XORNew(evalCopy, x0[i], x1[i])
			y1[i] = XORNew(evalCopy, x1[i], x2[i])
			y2[i] = XORNew(evalCopy, x2[i], x3[i])
			y3[i] = XORNew(evalCopy, x3[i], x0[i])
			ch <- true
		}(i)
	}
	for i := 0; i < 32; i++ {
		<-ch
	}

	for i := 0; i < 32; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			z0[i] = XORNew(evalCopy, y1[i], x3[i])
			z1[i] = XORNew(evalCopy, y2[i], x0[i])
			z2[i] = XORNew(evalCopy, y3[i], x1[i])
			z3[i] = XORNew(evalCopy, y0[i], x2[i])
			ch <- true
		}(i)
	}
	for i := 0; i < 32; i++ {
		<-ch
	}

	GF2FieldMul(aes.Evaluator, y0)
	GF2FieldMul(aes.Evaluator, y1)
	GF2FieldMul(aes.Evaluator, y2)
	GF2FieldMul(aes.Evaluator, y3)

	for i := 0; i < 32; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			z0[i] = XORNew(evalCopy, z0[i], y0[i])
			z1[i] = XORNew(evalCopy, z1[i], y1[i])
			z2[i] = XORNew(evalCopy, z2[i], y2[i])
			z3[i] = XORNew(evalCopy, z3[i], y3[i])
			ch <- true
		}(i)
	}
	for i := 0; i < 32; i++ {
		<-ch
	}
	z0 = append(z0, z1...)
	z0 = append(z0, z2...)
	z0 = append(z0, z3...)
	copy(x, z0)

}

func (aes *Ascon) ShiftRow(x []*rlwe.Ciphertext) {
	for i := 0; i < 8; i++ {
		x[1*8+i], x[5*8+i], x[9*8+i], x[13*8+i] = x[5*8+i], x[9*8+i], x[13*8+i], x[1*8+i]
		x[2*8+i], x[10*8+i], x[14*8+i], x[6*8+i] = x[10*8+i], x[14*8+i], x[6*8+i], x[2*8+i]
		x[3*8+i], x[15*8+i], x[11*8+i], x[7*8+i] = x[15*8+i], x[11*8+i], x[7*8+i], x[3*8+i]
	}
}

func (aes *Ascon) AddWhiteKey( pt, key []*rlwe.Ciphertext)  []*rlwe.Ciphertext {
	ch := make(chan bool, 128)
	for i := 0; i < 128; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			pt[i] = XORNew(evalCopy, pt[i], key[i])
			evalCopy.ScaleUp(pt[i], rlwe.NewScale(2), pt[i])
			ch <- true
		}(i)
	}
	for i := 0; i < 128; i++ {
		<-ch
	}	
	return pt
}

func (aes *Ascon) AddRoundKey( state, key []*rlwe.Ciphertext) {
	ch := make(chan bool, 128)
	for i := 0; i < 128; i++ {
		go func(i int) {
			evalCopy := aes.Evaluator.ShallowCopy()
			XOR(evalCopy, state[i], key[i], state[i])
			ch <- true
		}(i)
	}
	for i := 0; i < 128; i++ {
		<-ch
	}
	return
}

// Free XOR components**************************************
func (aes *Ascon) FreeAddRoundKey( state, key []*rlwe.Ciphertext) {
	for i := 0; i < len(state); i++ {
		FXOR(aes.Evaluator, state[i], key[i], state[i])
	}
}

func (aes *Ascon) FreeMixColumn( x []*rlwe.Ciphertext ) {
	x0, x1, x2, x3 := []*rlwe.Ciphertext{}, []*rlwe.Ciphertext{}, []*rlwe.Ciphertext{}, []*rlwe.Ciphertext{}
	
	for i := 0; i < 128; i++ {
		mod := i % 32
		if mod < 8 {
			x0 = append(x0, x[i])
		} else if mod >= 8 && mod < 16 {
			x1 = append(x1, x[i])
		} else if mod >= 16 && mod < 24 {
			x2 = append(x2, x[i])
		} else {
			x3 = append(x3, x[i])
		}
	}

	y0, y1, y2, y3 := make([]*rlwe.Ciphertext, len(x0)), make([]*rlwe.Ciphertext, len(x1)), make([]*rlwe.Ciphertext, len(x2)), make([]*rlwe.Ciphertext, len(x3))
	z0, z1, z2, z3 := make([]*rlwe.Ciphertext, len(x0)), make([]*rlwe.Ciphertext, len(x1)), make([]*rlwe.Ciphertext, len(x2)), make([]*rlwe.Ciphertext, len(x3))
	
	for i := 0; i < 32; i++ {
		y0[i] = FXORNew(aes.Evaluator, x0[i], x1[i])
		y1[i] = FXORNew(aes.Evaluator, x1[i], x2[i])
		y2[i] = FXORNew(aes.Evaluator, x2[i], x3[i])
		y3[i] = FXORNew(aes.Evaluator, x3[i], x0[i])
	}

	for i := 0; i < 32; i++ {
		z0[i] = FXORNew(aes.Evaluator, y1[i], x3[i])
		z1[i] = FXORNew(aes.Evaluator, y2[i], x0[i])
		z2[i] = FXORNew(aes.Evaluator, y3[i], x1[i])
		z3[i] = FXORNew(aes.Evaluator, y0[i], x2[i])
	}

	FreeGF2FieldMul(aes.Evaluator, y0)
	FreeGF2FieldMul(aes.Evaluator, y1)
	FreeGF2FieldMul(aes.Evaluator, y2)
	FreeGF2FieldMul(aes.Evaluator, y3)

	for i := 0; i < 32; i++ {
		z0[i] = FXORNew(aes.Evaluator, z0[i], y0[i])
		z1[i] = FXORNew(aes.Evaluator, z1[i], y1[i])
		z2[i] = FXORNew(aes.Evaluator, z2[i], y2[i])
		z3[i] = FXORNew(aes.Evaluator, z3[i], y3[i])
	}

	z0 = append(z0, z1...)
	z0 = append(z0, z2...)
	z0 = append(z0, z3...)
	copy(x, z0)
}