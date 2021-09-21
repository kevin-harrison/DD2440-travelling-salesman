import unittest
# Example input: (ᆨx1 ∨ x2 ∨  x3) ∧ (x2 ∨ x1 ∨ ᆨx3) = [(1,0,0, 1, 2, 3), (0,0,1, 2, 1 ,3)]


def calcExpectedValue(clauses, variables):
    expected_value = 0
    for (neg1, neg2, neg3, x1, x2, x3) in clauses: # Go through every clause encoded with form (ᆨ1, ᆨ2, ᆨ3, x1, x2, x3)

        # Calculate logical values of x1, x2, x3 when XORed with their negations. Set as None if they are still a random variable
        x1_val = None if variables[x1] == None else neg1 != variables[x1]
        x2_val = None if variables[x2] == None else neg2 != variables[x2]
        x3_val = None if variables[x3] == None else neg3 != variables[x3]

        # If any value in the clause is true then clause is worth 1 "point"
        if x1_val or x2_val or x3_val:
            expected_value += 1
            continue

        # Otherwise calculate expected value of clause
        # Example: (False, x4, x8)
        # ==> 3 possible combinations != (False, False, False) with 4 possible combinations 
        # ==> Expected value of clause = 3/4
        # DOESNT WORK WHEN WE HAVE FOR EXAMPLE: (False, not x4, x4)
        num_rand_variables = sum(1 for val in (x1_val, x2_val, x3_val) if val == None)
        expected_value += ((2**num_rand_variables) - 1) / (2**num_rand_variables)
        
    return expected_value


def Max3SAT(logical_statement, m, n):
    variables = [None] * n              # Initialize variables where None signifies a random variable
    for i in range(len(variables)):     # Greedy iteration over setting variables to True/False

        # Expected value when x_i = True
        variables[i] = 1
        expected_value_true = calcExpectedValue(logical_statement, variables)

        # Expected value when x_i = False
        variables[i] = 0
        expected_value_false = calcExpectedValue(logical_statement, variables)

        if expected_value_true > expected_value_false:
            variables[i] = 1
    
    # Expected value once all variables are set
    return calcExpectedValue(logical_statement, variables)



# Tests
class test(unittest.TestCase):
    def testSingleClause(self):
        test_input = [(0,0,0, 0,1,2)]
        test_vars = [1, None, None]
        self.assertEquals(calcExpectedValue(test_input, test_vars), 1)

    def testNegatedClause(self):
        test_input = [(1,1,0, 0,1,2)]
        test_vars = [1, None, None]
        self.assertEquals(calcExpectedValue(test_input, test_vars), 0.75)

    def testVarIsTrueAndFalse(self):
        test_input = [(1,0,1, 0,0,2)]
        test_vars = [1, None, None]
        self.assertEquals(calcExpectedValue(test_input, test_vars), 1)

    def testMultiClause(self):
        test_input = [(1,0,0, 2,0,2), (1,1,1, 0,1,2)]
        test_vars = [1, 1, None]
        self.assertEquals(calcExpectedValue(test_input, test_vars), 1.5)

    def testRandVarIsRepeated(self):
        test_input = [(0,0,0, 2,1,1)]
        test_vars = [1, None, None, None]
        self.assertEquals(calcExpectedValue(test_input, test_vars), 0.5)

    def testRandVarIsTrueAndFalse(self):
        test_input = [(0,1,0, 2,2,3)]
        test_vars = [1, None, None, None]
        self.assertEquals(calcExpectedValue(test_input, test_vars), 1)



if __name__ == '__main__':
    unittest.main()