#!/usr/bin/env python3
# Python matrix library

class Matrix:
    """ Rows and columns of the matrix """
    rows = 0
    cols = 0
    
    """ 2-dimensional list of entries in the matrix """
    entries = [[]]
    
    def __init__(self, rows, cols):
        """ Initialize zero matrix with given size """
        self.rows = rows
        self.cols = cols
        
        self.entries = [[0 for i in range(cols)] for i in range(rows)]
    
    def identity(n):
        """ Return the nxn identity matrix """
        i = Matrix(n, n)
        for row in range(i.rows):
            i.set(row, row, 1)
        
        return i
    
    def at(self, row, col):
        """ Return the value of the entry at the given position """
        return self.entries[row][col]
    
    def set(self, row, col, value):
        """ Set the entry at the given row and column to the given value """
        self.entries[row][col] = value
    
    def printMatrix(self):
        """ Display the matrix """
        for row in range(self.rows):
            for col in range(self.cols):
                print(self.at(row, col), end=' ')
            print()
    
    def copy(self):
        """ Return a Matrix that is the same as the current Matrix """
        newMatrix = Matrix(self.rows, self.cols)
        
        for row in range(self.rows):
            for col in range(self.cols):
                newMatrix.set(row, col, self.at(row, col))
        
        return newMatrix
    
    def subMatrix(self, delRow, delCol):
        """ Return this matrix but with the given row and column deleted """
        s = Matrix(self.rows - 1, self.cols - 1)
        
        row = 0
        sRow = 0
        while row < self.rows:
            if row == delRow:
                row += 1
            
            col = 0
            sCol = 0
            while col < self.cols:
                if col == delCol:
                    col += 1
                
                if (col < self.cols) and (row < self.rows):
                    s.set(sRow, sCol, self.at(row, col))
                
                col += 1
                sCol += 1
            
            row += 1
            sRow += 1
        
        return s
    
    def determinant(self):
        """ Compute the determinant of a square matrix """
        if self.rows != self.cols:
            raise Exception('Not a square matrix')
        
        # Trivial case: 1x1 matrix
        if (self.rows == 1) and (self.cols == 1):
            return self.at(0, 0)
        
        det = 0
        sign = 1 # Alternates for each column
        
        for col in range(self.cols):
            det += sign * self.at(0, col) * self.subMatrix(0, col).determinant()
            sign *= -1
        
        return det
    
    def rref(self, returnRank = False):
        """
        Perform Gauss-Jordan elimination to put this matrix in reduced row echelon form.
        Currently only works for matrices with full rank and whose top-leftmost
        square submatrix is invertible.
        """
        r = self.copy()
        #pivots = []
        pivot = 0
        
        if r.rows >= r.cols: maxRank = r.cols # Maximum possible rank of r
        else: maxRank = r.rows
        
        for row in range(r.rows):
            if pivot >= maxRank: break
            
            # Ensure that the current pivot does not contain a 0
            if r.at(row, pivot) == 0:
                # Try exchanging the current row with one below
                exchangeRow = row + 1
                while exchangeRow < r.rows:
                    if r.at(exchangeRow, pivot) != 0:
                        # Exchange current row with exchangeRow
                        p = Matrix.identity(r.rows) # Permutation matrix
                        
                        p.set(row, row, 0)
                        p.set(row, exchangeRow, 1)
                        p.set(exchangeRow, exchangeRow, 0)
                        p.set(exchangeRow, row, 1)
                        
                        r = p * r;
                        break
                    
                    exchangeRow += 1
            
            if r.at(row, pivot) != 0:
                #pivots.append(pivot)
                for opRow in range(r.rows):
                    if opRow == row: continue
                    
                    scalar = -1 * r.at(opRow, pivot) / r.at(row, pivot)
                    for col in range(r.cols):
                        r.set(opRow, col, r.at(opRow, col) + (scalar * r.at(row, col)))
                    
            pivot += 1
        
        rank = 0
        for row in range(r.rows):
            nonZero = False
            for col in range(r.cols):
                if r.at(row, col) != 0:
                    nonZero = True
                    
                    coef = r.at(row, col)
                    while col < r.cols:
                        r.set(row, col, r.at(row, col) / coef)
                        col += 1
                    break
            
            if nonZero: rank += 1
        
        return r
    
    def getRow(self, row):
        """ Return list of values in row of this matrix """
        rowVector = []
        for col in range(self.cols):
            rowVector.append(self.at(row, col));
        return rowVector
    
    def getCol(self, col):
        """ Return list of values in column of this matrix """
        colVector = []
        for row in range(self.rows):
            colVector.append(self.at(row, col))
        return colVector
    
    def scalarMult(self, c):
        """ Return this matrix multiplied by a scalar c """
        ca = Matrix(self.rows, self.cols)
        for row in range(self.rows):
            for col in range(self.cols):
                ca.set(row, col, c * self.at(row, col))
        
        return ca
    
    def __add__(a, b):
        """ Return the sum of a and b """
        if (a.cols != b.cols) or (a.rows != b.rows):
            raise Exception('Matrices are of incompatible sizes')
        
        c = Matrix(a.rows, a.cols)
        for row in range(a.rows):
            for col in range(a.cols):
                c.set(row, col, a.at(row, col) + b.at(row, col))
        
        return c
    
    def __sub__(a, b):
        """ Return the difference of a and b """
        return a + (b.scalarMult(-1))
    
    def __mul__(a, b):
        """ Return the matrix product of a and b """
        if a.cols != b.rows:
            raise Exception('Matrices are of incompatible sizes')
        
        ab = Matrix(a.rows, b.cols)
        
        for row in range(ab.rows):
            for col in range(ab.cols):
                aRow = a.getRow(row)
                bCol = b.getCol(col)
                value = sum([i*j for i, j in zip(aRow, bCol)])
                ab.set(row, col, value)
        
        return ab

# Run simple tests/sanity checks
if __name__ == '__main__':
    values = [[1, 1, 2],
              [2, 2, 1]]
    
    a = Matrix(len(values), len(values[0]))
    i = 0
    for row in values:
        j = 0
        for col in row:
            a.set(i, j, col)
            j += 1
        i += 1
    
    a.printMatrix()
    print('------------')
    
    a.rref().printMatrix()
    print('------------')
    
    