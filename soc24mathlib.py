"""ASSIGNMENT 1- """

"""Function 1. Returns the GCD of a and b"""
def pair_gcd(a: int, b: int) -> int:   #Based on Euclidian algorithm
    if not a or not b:
        raise ValueError("Inputs must be non-zero integers.")
    if a < b:
        a, b = b, a
    rem = a % b
    while (rem != 0):
        a = b
        b = rem
        rem = a % b
    return b


def pair_egcd(a: int, b: int) -> tuple[int, int, int]:
    """
    This function implements the Extended Euclidean Algorithm recursively to compute
    the greatest common divisor (GCD) of two integers a and b, along with
    the Bezout coefficients x and y such that ax + by = gcd(a, b).

    Args:
        a: The first integer.
        b: The second integer.

    Returns:
        A tuple (x, y, d) where d is the GCD of a and b, and x and y are
        integers such that ax + by = d.
    """
    if a == 0 and b == 0:
        raise ValueError("Both a and b cannot be zero")

    if b == 0:
        return (1, 0, a)

    x1, y1, d = pair_egcd(b, a % b)
    x, y = y1, x1 - (a // b) * y1
    return (x, y, d)

"""Function 3: Returns the GCD of all integers provided as arguments"""
def gcd(*args: int) -> int:
  """
  This function implements the GCD of all integers provided as arguments using the pair_gcd function.

  Args:
      *args: Variable number of integers.

  Returns:
      The greatest common divisor (GCD) of all the input integers.
  """

  if not args:
    raise ValueError("No arguments provided")

  # Initialize with the first argument
  result = args[0]
  for num in args[1:]:
    # Use pair_gcd to find GCD between current result and each argument
    result = pair_gcd(result, num)

  return result


"""Function 4: Returns the LCM of a and b"""
def pair_lcm(a: int, b: int) -> int:
  """
  This function implements the Least Common Multiple (LCM) of two integers a and b.

  Args:
      a: The first integer.
      b: The second integer.

  Returns:
      The least common multiple (LCM) of a and b.
  """

  if not a or not b:
    raise ValueError("Both a and b cannot be zero")

  gcd = pair_gcd(a, b)
  # Use the recursive property of LCM: lcm(a, b) = (a * b) // gcd(a, b)
  return a * b // gcd


"""Function 5: Returns the LCM of all integers provided as arguments"""
def lcm(*args: int) -> int:
  """
  This function implements the Least Common Multiple (LCM) of all integers provided as arguments using the pair_lcm function.

  Args:
      *args: Variable number of integers.

  Returns:
      The least common multiple (LCM) of all the input integers.
  """

  if not args:
    raise ValueError("No arguments provided")

  # Initialize with the first argument
  result = args[0]
  for num in args[1:]:
    # Use pair_lcm to find LCM between current result and each argument
    result = pair_lcm(result, num)

  return result


"""Function 6: Returns True if a and b are relatively prime, False otherwise"""
def are_relatively_prime(a: int, b: int) -> bool:
  """
  This function checks if two integers a and b are relatively prime (coprime) using the GCD property.

  Args:
      a: The first integer.
      b: The second integer.

  Returns:
      True if a and b are relatively prime (coprime), False otherwise.
  """
  if not a or not b:
    raise ValueError("Both a and b cannot be zero")
  return pair_gcd(a, b) == 1  # GCD of coprime numbers is always 1


def mod_inv(a: int, n: int) -> int:
    """
    This function computes the modular inverse of a modulo n using the Extended Euclidean Algorithm (EEA).

    Args:
        a: The integer for which to find the modular inverse.
        n: The modulus.

    Returns:
        The modular inverse of a modulo n, or raises a ValueError if the inverse doesn't exist.
    """
    if not a or not n:
        raise ValueError("Both a and n cannot be zero")

    x, y, d = pair_egcd(a, n)
    
    # Debug statements
    # print(f"pair_egcd({a}, {n}) = (x={x}, y={y}, d={d})")

    # Check if modular inverse exists (GCD of a and n must be 1)
    if d != 1:
        raise ValueError("Modular inverse does not exist (a and n are not coprime)")

    # Modular inverse correction
    inv = x % n
    # print(f"Modular inverse of {a} mod {n} is {inv}")
    return inv




"""Function 8: Returns the unique value of a modulo product of all n[i] such that a = a[i] (mod n[i]) using Chinese Remainder Theorem"""
def crt(a: list[int], n: list[int]) -> int:
  """
  This function implements the Chinese Remainder Theorem (CRT) to find the unique solution 
  of the system of congruences a = a[i] (mod n[i]) for all i.

  Args:
      a: A list of integers representing residues.
      n: A list of integers representing moduli (all pairwise coprime).

  Returns:
      The unique solution of the system of congruences.

  Raises:
      ValueError: If the lengths of a and n don't match or if the moduli are not pairwise coprime.
  """

  if len(a) != len(n):
    raise ValueError("Length of a and n must be equal")

  # Check for pairwise coprimality using pair_gcd
  for i in range(len(n)):
    for j in range(i + 1, len(n)):
      if pair_gcd(n[i], n[j]) != 1:
        raise ValueError("Moduli must be pairwise coprime")

  # Pre-calculate products and inverses
  products = [1] * len(n)
  inverses = [0] * len(n)
  for i in range(len(n)):
    products[i] = prod(n[:i] + n[i + 1:])  # Product of all other moduli
    inverses[i] = mod_inv(products[i], n[i])  # Modular inverse of product

  # Calculate the solution using CRT formula
  result = 0
  for i in range(len(n)):
    result = (result + a[i] * products[i] * inverses[i]) % (prod(n))
  return result

# Helper function to calculate product
def prod(iterable):
  product = 1
  for num in iterable:
    product *= num
  return product



"""Function 9: Returns a^m (mod n) computed using a method called fast exponentiation."""
def pow(a: int, m: int, n: int) -> int:
  """
  This function implements fast exponentiation to calculate a^m modulo n.

  Args:
      a: The base of the exponentiation.
      m: The exponent.
      n: The modulus.

  Returns:
      The value of a^m modulo n.
  """

  if not m:
    return 1

  result = 1
  base = a % n  

  while m > 0:
    if m & 1:
      result = (result * base) % n  # Include base if LSB of m is 1
    m >>= 1  # Right shift m by 1
    base = (base * base) % n  # Square base modulo n

  return result

def is_quadratic_residue_prime(a: int, p: int) -> int:
    """
    Determines the quadratic residue status of an integer a modulo p using Euler's Criterion.

    Args:
        a (int): The integer to be evaluated.
        p (int): The prime modulus.

    Returns:
        1: If a is a quadratic residue modulo p.
        -1: If a is a quadratic non-residue modulo p.
        0: If a is not coprime to p.
    """
    if not (a and p):
        return 0
    
    # Check if a is coprime to p
    if not are_relatively_prime(a, p):
        return 0

    # Reduce a modulo p for efficiency
    a %= p

    # Euler's Criterion
    result = pow(a, (p - 1) // 2, p)
    return 1 if result == 1 else -1

def is_quadratic_residue_prime_power(a: int, p: int, e: int) -> int:
    """
    Determines the quadratic residue status of an integer a modulo p^e using properties of quadratic residues.

    Args:
        a (int): The integer to be checked.
        p (int): The base prime number.
        e (int): The exponent to which the prime is raised.

    Returns:
        1: If a is a quadratic residue modulo p^e.
        -1: If a is a quadratic non-residue modulo p^e.
        0: If a is not coprime to p^e.
    """
    if not (a and p):
        return 0

    # Check if a is coprime to p^e
    if not are_relatively_prime(a, p**e):
        return 0

    # Reduce a modulo p^e for efficiency
    a %= p**e

    # Special case for e = 1 (reduce to standard check)
    if e == 1:
        return is_quadratic_residue_prime(a, p)

    # Use properties of quadratic residues
    residue_check = is_quadratic_residue_prime(a, p)
    if residue_check == -1:
        return -1

    # If a is non-residue mod p and e is even, it's a residue mod p^e
    return 1 if e % 2 == 0 else -1
