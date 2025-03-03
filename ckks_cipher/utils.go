package ckks_cipher

import (
	"fmt"
	"math"
	"strings"

	"crypto/rand"
	"io"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"golang.org/x/crypto/sha3"
)

type BitSet struct {
    bits []uint8
    size  int
}


func NewBitSet(size int) *BitSet {
    return &BitSet{
        bits: make([]uint8, size),
        size:  size,
    }
}

func (b *BitSet) Set(x int) {
	if x<0 {panic("x value must > 0")}
	for i:=0;i<b.size;i++{
		b.bits[i] = uint8(x & 0x1)
		x >>= 1
	}
}

func (b *BitSet) ToULong() uint64 {
	var Out uint64
	for i, bit := range b.bits {
		Out = Out + uint64(bit << i)
	}
	return Out
}

func (b *BitSet) ToString() string {
    var sb strings.Builder
    for _, byte := range b.bits {
        if byte == 1 {
            sb.WriteString("1")
        } else {
            sb.WriteString("0")
        }
    }
    return sb.String()
}

func (b *BitSet) Copy() *BitSet {
    out := NewBitSet(b.size)
    for i, bit := range b.bits {
        if bit == 1 {
            out.bits[i] = 1
        } else {
            out.bits[i] = 0
        }
    }
    return out
}

func Xor(v0, v1 *BitSet) *BitSet {
	
    if v0.size!= v1.size {
        panic("bit sets have different sizes")
    }
	Out := NewBitSet(v0.size)
    for i := 0; i < v0.size; i++ {
        Out.bits[i] = v0.bits[i] ^ v1.bits[i]
    }
	return Out
}

func (b *BitSet) GetSize() int {
    return b.size
}

func vectorLeftRotationInplace(vector []*rlwe.Ciphertext, rotNum int) {
    for i := 0; i < rotNum; i++ {
        first := vector[0]
        copy(vector, vector[1:])
        vector[len(vector)-1] = first
    }
}

func vectorRightRotationInplace(vector []*rlwe.Ciphertext, rotNum int) {
    for i := 0; i < rotNum; i++ {
        last := vector[len(vector)-1]
        copy(vector[1:], vector)
        vector[0] = last
    }
}

func ctr(iv *BitSet, ctr uint64) *BitSet {
	out := iv.Copy()
	for i := iv.size - 64; i < iv.size; i++ {
		out.bits[i] = iv.bits[i] ^ byte( (ctr>>uint(i-8))&1 )
	}
	return out
}
// MinInt returns the minimum value of the input of int values.
func MinInt(a, b int) (r int) {
	if a <= b {
		return a
	}
	return b
}


func rankOfMatrix(matrix [][]bool) uint {
	mat := make([][]bool, len(matrix))
	copy(mat, matrix)

	size := uint(len(mat[0]))
	var row uint = 0

	// Transform to upper triangular matrix
	for col := uint(1); col <= size; col++ {
		if !mat[row][size-col] {
			r := row
			for r < uint(len(mat)) && !mat[r][size-col] {
				r++
			}
			if r >= uint(len(mat)) {
				continue
			} else {
				mat[row], mat[r] = mat[r], mat[row]
			}
		}
		for i := row + 1; i < uint(len(mat)); i++ {
			if mat[i][size-col] {
				for j := uint(0); j < size; j++ {
					mat[i][j] = mat[i][j] != mat[row][j] // XOR operation
				}
			}
		}
		row++
		if row == size {
			break
		}
	}
	return row
}

func generateRandomMatrix(rows, cols int) [][]bool {
	// Calculate the required bytes, ensuring to round up for any extra bits
	requiredBytes := (rows*cols + 7) / 8 // Round up to nearest byte

	// Initialize SHAKE256 with a secure random seed
	hash := make([]byte, 64) // Use a 512-bit seed
	io.ReadFull(rand.Reader, hash)
	shake := sha3.NewShake256()
	shake.Write(hash)

	// Generate enough random bytes
	randomBytes := make([]byte, requiredBytes)
	shake.Read(randomBytes)

	matrix := make([][]bool, rows)
	for i := range matrix {
		matrix[i] = make([]bool, cols)
		for j := 0; j < cols; j++ {
			byteIndex := (i*cols + j) / 8
			bitIndex := uint((i*cols + j) % 8)
			matrix[i][j] = (randomBytes[byteIndex] & (1 << bitIndex)) != 0
		}
	}
	return matrix
}

// generateRandomVector creates a random vector of booleans with the provided length using SHAKE256
func generateRandomVector(length int) []bool {
	// Calculate the number of bytes needed to hold `length` bits
	requiredBytes := (length + 7) / 8 // Round up to ensure entire vector fits

	// Initialize SHAKE256 with a secure random seed
	hash := make([]byte, 64) // Use a 512-bit seed
	io.ReadFull(rand.Reader, hash)
	shake := sha3.NewShake256()
	shake.Write(hash)

	// Generate enough random bytes
	randomBytes := make([]byte, requiredBytes)
	shake.Read(randomBytes)

	vector := make([]bool, length)
	for i := 0; i < length; i++ {
		byteIndex := i / 8
		bitIndex := uint(i % 8)
		vector[i] = (randomBytes[byteIndex] & (1 << bitIndex)) != 0
	}
	return vector
}

// ToFloatMat converts a matrix of boolean values to a matrix of float64 values
func ToFloatMat(boolMat [][]bool) [][]float64 {
	// Create a new matrix of the same size to hold float64 values
	matrix := make([][]float64, len(boolMat))
	
	for i := 0; i < len(boolMat); i++ {
		// Initialize each row of the float64 matrix
		matrix[i] = make([]float64, len(boolMat[i]))
		
		for j := 0; j < len(boolMat[i]); j++ {
			// Convert each boolean value to float64 (true -> 1.0, false -> 0.0)
			if boolMat[i][j] {
				matrix[i][j] = 1.0
			} else {
				matrix[i][j] = 0.0
			}
		}
	}
	
	return matrix
}

// ToFloatVec converts a vector of boolean values to a vector of float64 values
func ToFloatVec(boolVec []bool) []float64 {
	// Create a new slice to hold float64 values with the same length as the boolean vector
	floatVec := make([]float64, len(boolVec))

	// Convert each boolean value to float64 (true -> 1.0, false -> 0.0)
	for i, val := range boolVec {
		if val {
			floatVec[i] = 1.0
		} else {
			floatVec[i] = 0.0
		}
	}

	return floatVec
}


// func main() {
// 	rows, cols := 351, 351 // Example size
// 	matrix := generateRandomMatrix(rows, cols)

// 	// Display the generated matrix
// 	for _, row := range matrix {
// 		for _, val := range row {
// 			if val {
// 				fmt.Print("1 ")
// 			} else {
// 				fmt.Print("0 ")
// 			}
// 		}
// 		fmt.Println()
// 	}

// 	rank := rankOfMatrix(matrix)
//     vec := generateRandomVector(351)
//     fmt.Println(vec)
// 	fmt.Printf("Rank of the matrix: %d\n", rank)
// }

// PrintPrecisionStats decrypts, decodes and prints the precision stats of a ciphertext.
func PrintPrecisionStatsBoolean(params ckks.Parameters, ct *rlwe.Ciphertext, ecd *ckks.Encoder, dec *rlwe.Decryptor) {

	var err error

	// Decrypts the vector of plaintext values
	pt := dec.DecryptNew(ct)

	// Decodes the plaintext
	have := make([]float64, params.MaxSlots())
	if err = ecd.Decode(pt, have); err != nil {
		panic(err)
	}
    want := make([]float64, params.MaxSlots())

	for i:=0;i<len(have);i++{
		value := math.Abs(have[i])
		if ( math.Abs( value - want[i] ) > 0.1){
            want[i] = 1.0
		} else {
            want[i] = 0.0
        }
	}

	// Pretty prints some values
	fmt.Printf("Have: ")
	for i := 0; i < 4; i++ {
		fmt.Printf("%20.15f ", have[i])
	}
	fmt.Printf("...\n")

	fmt.Printf("Want: ")
	for i := 0; i < 4; i++ {
		fmt.Printf("%20.15f ", want[i])
	}
	fmt.Printf("...\n")

	// Pretty prints the precision stats
	fmt.Println(ckks.GetPrecisionStats(params, ecd, dec, have, want, 0, false).String())
}
