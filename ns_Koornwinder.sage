import itertools
import numpy
from sage.structure.list_clone import ClonableElement, ClonableList, ClonableIntArray
from sage.structure.parent import Parent

def alpha_check(i, n):
    r"""
    returns the ith simple coroot in the CC_n system.
    """
    if i == 0:
        return( tuple([-1] + (n-1)*[0] + [1/2]) )
    elif i == n:
        return( tuple((n-1)*[0] + [1, 0]) )
    else:
        return( tuple((i-1)*[0] + [1, -1] + (n-i)*[0]) )

def coroot(a, b, K, n):
    r"""
    returns type CC_n coroot with finite component determined by pair (a, b)
    """
    if a == 0 or abs(a) > n or abs(b) > n:
        raise(ValueError("root coordinates out of bounds"))
    try:
        lst_form = [Integer(a / abs(a)) if i == abs(a) else Integer(b/abs(b)) if i == abs(b) else 0 for i in range(1, n+1)]
    except:
        raise(TypeError(f"the numbers are f{(a, b, K)}, with type={(type(v) for v in [a, b, K])}"))
    lst_form += [Integer(2*K)/2]
    return(tuple(lst_form))
    

class KoornwinderBoxDiagram(ClonableList):
    r"""
    Class for Koornwinder box diagrams
    """

    def __init__(self, mu):
        ClonableList.__init__(self, parent = Parent(), lst = mu)
        self.check()

    def check(self):
        r"""checks whether mu is a valid CC weight"""
        if not([type(v) == int for v in self]):
            raise ValueError("Values are not integers")

    def _compute_v_mu(self):
        r"""
        Initialize the signed permutation v_mu as a parameter.  Access by self.v_mu().
        """
        v_mu = self.__len__() * [1]
        for i in range(self.__len__()):
            i_abs = abs(self.__getitem__(i))
            for j in range(i+1, self.__len__()):
                j_abs = abs(self.__getitem__(j))
                if i_abs < j_abs:
                    v_mu[i] += 1
                elif i_abs > j_abs:
                    v_mu[j] += 1
                elif i_abs == j_abs:
                    if self.__getitem__(i) > 0:
                        v_mu[i] += 1
                    elif self.__getitem__(i) <= 0:
                        v_mu[j] += 1
        for i in range(self.__len__()):
            if self.__getitem__(i) > 0:
                v_mu[i] *= -1
        self._v_mu = v_mu
    
    def v_mu(self):
        r"""
        Retrieval method for the signed permutation v_mu which takes self into antidominant chamber

        EXAMPLES:: 
        sage: D = KoornwinderBoxDiagram([-3, 5, -1, 4])
        sage: D.v_mu()
        [3, -1, 4, -2]

        sage: D = KoornwinderBoxDiagram([0, 2, 3, -1, 1])
        sage: D.v_mu()
        [5, -2, -1, 3, -4] 
        """
        try:
            return(self._v_mu)
        except:
            self._compute_v_mu()
            return(self._v_mu)
    
    def dg_plus(self):
        r"""
        alias of self.boxes()
        """

    def boxes(self):
        r"""
        Returns a list of boxes in self, represented as coordinate pairs (r, s), 
        and in the spiral reading order.

        EXAMPLES::  
        sage: D = KoornwinderBoxDiagram([-3, 5, -1, 4])
        sage: D.pretty_print_boxes()
         _ _ _            
        |_|_|_|0|_ _ _ _ _  
             _|1|_|_|_|_|_| 
            |_|2|_ _ _ _    
              |3|_|_|_|_|   
        sage: D.boxes()
        [(1, 1), (3, 1), (2, -1), (0, -1), (1, 2), (3, 2), (0, -2), (1, 3), (3, 3), (0, -3), 
        (1, 4), (3, 4), (1, 5)]

        sage: D = KoornwinderBoxDiagram([0, 2, 3, -1, 1])
        sage: D.pretty_print_boxes()
          |0|_ _    
          |1|_|_|_  
         _|2|_|_|_| 
        |_|3|_      
          |4|_| 
        sage: D.boxes()
        [(1, 1), (2, 1), (4, 1), (3, -1), (1, 2), (2, 2), (2, 3)]
        """
        out = []
        for r in range(len(self)):
            C = self[r]
            if C >= 0:
                for c in range(C):
                    out.append((r, c+1))
            if C < 0:
                for c in range(1, 1-C):
                    out.append((r, -c))
        out.sort(key = (lambda x : (x[1], x[0]) if x[1] > 0 else (-x[1], 2*self.__len__() - x[0])))
        return(out)

    def interior_box(self, r, c):
        r"""
        Returns True if self has a box at position (r, c) and False otherwise.  
        Defaults True when c == 0.

        EXAMPLES::  
        sage: D = KoornwinderBoxDiagram([-3, 5, -1, 4])
        sage: D.pretty_print_boxes()
         _ _ _            
        |_|_|_|0|_ _ _ _ _  
             _|1|_|_|_|_|_| 
            |_|2|_ _ _ _    
              |3|_|_|_|_|   
        sage: D.boxes()
        [(1, 1), (3, 1), (2, -1), (0, -1), (1, 2), (3, 2), (0, -2), (1, 3), (3, 3), (0, -3), 
        (1, 4), (3, 4), (1, 5)]
        sage: D.interior_box(2, 1)
        False
        sage: D.interior_box(3, 1)
        True
        sage: D.interior_box(0, 0)
        True

        sage: D = KoornwinderBoxDiagram([0, 2, 3, -1, 1])
        sage: D.pretty_print_boxes()
          |0|_ _    
          |1|_|_|_  
         _|2|_|_|_| 
        |_|3|_      
          |4|_| 
        sage: D.boxes()
        [(1, 1), (2, 1), (4, 1), (3, -1), (1, 2), (2, 2), (2, 3)]
        sage: D.interior_box(2, 1)
        True
        sage: D.interior_box(3, 1)
        False
        sage: D.interior_box(0, 0)
        True
        """
        r_val = self.__getitem__(r)
        if r_val * c < 0:
            return(False)
        elif abs(c) > abs(r_val):
            return(False)
        else:
            return(True)

    def _check_box_call(self, r, c):
        r"""
        Raises ValueError if (r, c) is not an interior box
        """
        if not self.interior_box(r, c):
            raise(ValueError(f"expected coordinates ({r}, {c}) to be within the shape mu"))

    def pretty_print_boxes(self):
        r"""
        Pretty print the argument in an ascii_art style.

        EXAMPLES:: 
        sage: D = KoornwinderBoxDiagram([-3, 5, -1, 4])
        sage: D.pretty_print_boxes()
         _ _ _            
        |_|_|_|0|_ _ _ _ _  
             _|1|_|_|_|_|_| 
            |_|2|_ _ _ _    
              |3|_|_|_|_|   

        sage: D = KoornwinderBoxDiagram([0, 2, 3, -1, 1])
        sage: D.pretty_print_boxes()
          |0|_ _    
          |1|_|_|_  
         _|2|_|_|_| 
        |_|3|_      
          |4|_|    
        """
        lst = self._get_list()
        transition_intervals = []
        bar_intervals = []
        lm = 0
        rm = 0
        for i in range(self.__len__()+1):
            if i == self.__len__():
                transition_intervals.append(sorted([0, self.__getitem__(-1)]))
                bar_intervals.append(sorted([0, self.__getitem__(-1)]))
            else:
                lm = min(self.__getitem__(i), lm)
                rm = max(self.__getitem__(i), rm)
                if i == 0:
                    transition_intervals.append(sorted([0, self.__getitem__(i)]))
                    bar_intervals.append([0, -1])
                else:
                    transition_intervals.append(sorted([self.__getitem__(i-1), self.__getitem__(i)]))
                    bar_intervals.append(sorted([0, self.__getitem__(i-1)]))
        print_string = ""
        for j in range(lm, rm+1):
            ti = transition_intervals[0]
            if ti[0] <= j and j < ti[1]:
                print_string += " _"
            else:
                print_string += "  "
        print_string += "\n"
        for i in range(self.__len__()):
            ti = transition_intervals[i+1]
            bi = bar_intervals[i+1]
            for j in range(lm, rm+1):
                if j == 0:
                    print_string += f"|{i}"
                if bi[0] <=j and j < bi[1]:
                    print_string += "|_"
                elif j == bi[1] and j < ti[1]:
                    print_string += "|_"
                elif j == bi[1] and j == ti[1]:
                    print_string += "| "
                elif ti[0] <= j and j < ti[1]:
                    print_string += " _"
                else:
                    print_string += "  "
            print_string += "\n"
        print(print_string)


    def attacking_boxes(self, r, c):
        """
        Returns a list of the boxes attacking (r, c).

        EXAMPLES:: 
        sage: D = KoornwinderBoxDiagram([-3, 5, -1, 4])
        sage: D.pretty_print_boxes()
         _ _ _            
        |_|_|_|0|_ _ _ _ _  
             _|1|_|_|_|_|_| 
            |_|2|_ _ _ _    
              |3|_|_|_|_|   
        sage: D.attacking_boxes(1, 2)
        [(0, -1), (2, -1), (3, 1)]
        sage: D.attacking_boxes(1, 3)
        [(0, -2), (3, 2)]

        sage: D = KoornwinderBoxDiagram([0, 2, 3, -1, 1])
        sage: D.pretty_print_boxes()
          |0|_ _    
          |1|_|_|_  
         _|2|_|_|_| 
        |_|3|_      
          |4|_|    

        sage: D.attacking_boxes(2, 1)
        [(1, 1)]
        sage: D.attacking_boxes(3, -1)
        [(4, 1), (2, 1), (1, 1)]
        """
        self._check_box_call(r, c)
        r_val = self.__getitem__(r)
        same_col = []
        opp_side = []
        prev_col = []
        if c > 0:
            for i in range(self.__len__()):
                i_val = self.__getitem__(i)
                if i_val == 0:
                    if c == 1 and i < r:
                        opp_side.append( (i, 0) )
                    else:
                        continue
                elif i_val < 0:
                    if c == 1:
                        if i < r:
                            opp_side.append( (i, 0) )
                    elif abs(c) - 1 <= abs(i_val):
                        opp_side.append( (i, 1-c) )
                elif c != 1 and i > r:
                    if abs(c) - 1 <= abs(i_val):
                        prev_col.append( (i, c-1) )
                elif i < r:
                    if abs(c) <= abs(i_val):
                        same_col.append( (i, c) )
            attack = list(reversed(same_col)) + opp_side + list(reversed(prev_col))
        elif c < 0:
            for i in range(self.__len__()):
                i_val = self.__getitem__(i)
                if i_val == 0:
                    if c == -1:
                        opp_side.append( (i, 1+c) )
                    else:
                        continue
                elif i_val > 0:
                    if abs(c) - 1 <= abs(i_val):
                        opp_side.append( (i, -c) )
                elif i < r:
                    if abs(c) - 1 <= abs(i_val):
                        prev_col.append( (i, c+1) )
                elif i > r:
                    if abs(c) <= abs(i_val):
                        same_col.append( (i, c) )
            attack = same_col + list(reversed(opp_side)) + prev_col
        return(attack)
    
    def leg_length(self, r, c):
        r"""
        returns the leg length statistic of the box (r, c)

        EXAMPLES:: 
        sage: D = KoornwinderBoxDiagram([-3, 5, -1, 4])
        sage: D.pretty_print_boxes()
         _ _ _            
        |_|_|_|0|_ _ _ _ _  
             _|1|_|_|_|_|_| 
            |_|2|_ _ _ _    
              |3|_|_|_|_|   
        sage: D.leg_length(1, 2)
        3
        sage: D.leg_length(1, 3)
        2

        sage: D = KoornwinderBoxDiagram([0, 2, 3, -1, 1])
        sage: D.pretty_print_boxes()
          |0|_ _    
          |1|_|_|_  
         _|2|_|_|_| 
        |_|3|_      
          |4|_|    

        sage: D.leg_length(2, 1)
        2
        sage: D.leg_length(3, -1)
        0
        """
        self._check_box_call(r, c)
        return(abs(self.__getitem__(r)) - abs(c))

    def bgrw_filling(self, r, c):
        r"""
        returns the portion of the box greedy reduced word coming from box (r, c)

        EXAMPLES:: 
        sage: D = KoornwinderBoxDiagram([-3, 5, -1, 4])
        sage: D.pretty_print_boxes()
         _ _ _            
        |_|_|_|0|_ _ _ _ _  
             _|1|_|_|_|_|_| 
            |_|2|_ _ _ _    
              |3|_|_|_|_|   
        sage: D.bgrw_filling(1, 2)
        [4, 3, 2, 1, 0]
        sage: D.bgrw_filling(1, 3)
        [3, 4, 3, 2, 1, 0]

        sage: D = KoornwinderBoxDiagram([0, 2, 3, -1, 1])
        sage: D.pretty_print_boxes()
          |0|_ _    
          |1|_|_|_  
         _|2|_|_|_| 
        |_|3|_      
          |4|_|    

        sage: D.bgrw_filling(2, 1)
        [2, 1, 0]
        sage: D.bgrw_filling(3, -1)
        [5, 4, 3, 2, 1, 0]
        """
        self._check_box_call(r, c)
        if abs(c) == 1:
            if c > 0:
                return([s for s in range(r, -1, -1)])
            elif c < 0:
                return([s for s in range(self.__len__(), -1, -1)])
        else:
            atk = len(self.attacking_boxes(r, c))
            return([s for s in range(atk+1, self.__len__())] + [s for s in range(self.__len__(), -1, -1)])
            

    def box_greedy_reduced_word(self):
        r"""
        returns box greedy reduced word for self as list

        EXAMPLES:: 
        sage: D = KoornwinderBoxDiagram([-3, 5, -1, 4])
        sage: D.pretty_print_boxes()
         _ _ _            
        |_|_|_|0|_ _ _ _ _  
             _|1|_|_|_|_|_| 
            |_|2|_ _ _ _    
              |3|_|_|_|_|   
        sage: D.box_greedy_reduced_word()
        [1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 4, 3, 2, 1, 0, 4, 3, 2, 1, 0, 4, 3, 2, 1, 0, 3, 4, 3, 
        2, 1, 0, 3, 4, 3, 2, 1, 0, 3, 4, 3, 2, 1, 0, 3, 4, 3, 2, 1, 0, 3, 4, 3, 2, 1, 0, 3, 4, 3, 
        2, 1, 0, 2, 3, 4, 3, 2, 1, 0]

        sage: D = KoornwinderBoxDiagram([0, 2, 3, -1, 1])
        sage: D.pretty_print_boxes()
          |0|_ _    
          |1|_|_|_  
         _|2|_|_|_| 
        |_|3|_      
          |4|_|    

        sage: D.box_greedy_reduced_word()
        [1, 0, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 4, 5, 4, 3, 2, 1, 0, 4, 5, 4, 3, 2, 1, 0, 
        1, 2, 3, 4, 5, 4, 3, 2, 1, 0]
        """
        return(list(itertools.chain.from_iterable([self.bgrw_filling(r, c) for (r, c) in self.boxes()])))

    def root_sequence(self, r, c):
        r"""
        returns root sequence of box (r, c)

        EXAMPLES:: 
        sage: D = KoornwinderBoxDiagram([0, 2, 3, -1, 1])
        sage: D.pretty_print_boxes()
          |0|_ _    
          |1|_|_|_  
         _|2|_|_|_| 
        |_|3|_      
          |4|_|    

        sage: D.root_sequence(1, 1)
        [(0, 1, 0, 0, 1, -2), (0, 1, 0, 0, 0, -3/2)]
        sage: D.root_sequence(2, 1)
        [(1, 0, 0, 0, 1, -3), (1, 1, 0, 0, 0, -4), (1, 0, 0, 0, 0, -5/2)]
        sage: D.root_sequence(4, 1)
        [(0, 0, 1, 1, 0, -2), (0, 0, 0, 1, 1, -1), (0, 1, 0, 1, 0, -2), (1, 0, 0, 1, 0, -3)
        (0, 0, 0, 1, 0, -1/2)]
        sage: D.root_sequence(3, -1)
        [(0, 0, 1, 0, 0, -1), (0, 0, 1, 0, 1, -1), (0, 1, 1, 0, 0, -2), (1, 0, 1, 0, 0, -3), 
        (0, 0, 1, 1, 0, -1), (0, 0, 1, 0, 0, -1/2)]
        sage: D.root_sequence(1, 2)
        [(0, 1, 0, 0, -1, -1), (0, 1, 0, 0, 0, -1), (0, 1, 0, 0, 1, -1), (1, 1, 0, 0, 0, -3), 
        (0, 1, 0, 1, 0, -1), (0, 1, 1, 0, 0, -1), (0, 1, 0, 0, 0, -1/2)]
        sage: D.root_sequence(2, 2)
        [(1, 0, 0, 0, -1, -2), (1, 0, 0, 0, 0, -2), (1, 0, 0, 0, 1, -2), (1, 0, 0, 1, 0, -2)
        (1, 0, 1, 0, 0, -2), (1, 1, 0, 0, 0, -2), (1, 0, 0, 0, 0, -3/2)]
        sage: D.root_sequence(2, 3)
        [(1, -1, 0, 0, 0, -1), (1, 0, -1, 0, 0, -1), (1, 0, 0, -1, 0, -1), (1, 0, 0, 0, -1, -1)
        (1, 0, 0, 0, 0, -1), (1, 0, 0, 0, 1, -1), (1, 0, 0, 1, 0, -1), (1, 0, 1, 0, 0, -1)
        (1, 1, 0, 0, 0, -1), (1, 0, 0, 0, 0, -1/2)]
        """
        self._check_box_call(r, c)
        start = -abs(self.v_mu()[r])
        atks = self.attacking_boxes(r, c)
        leg = self.leg_length(r, c)
        r_seq = []
        v_mu = self.v_mu()
        if c == 1:
            for ell in range(r-1, -1, -1):
                aj_r = atks[ell][0]
                aj_c = atks[ell][1]
                if aj_c > 0:
                    r_seq.append( (start, -abs(v_mu[aj_r]), (leg + self.leg_length(aj_r, aj_c) + 1)))
                else:
                    r_seq.append( (start, -abs(v_mu[aj_r]), (leg + self.leg_length(aj_r, aj_c) + 1)))
            #print(leg + 1/2)
            r_seq.append( (start, 0, (leg + 1/2)) )
            r_seq = [coroot(a, b, K, self.__len__()) for (a, b, K) in r_seq]
            return(r_seq)
        else:
            for ell in range(len(atks)+2, self.__len__()+1):
                r_seq.append( (start, ell, (leg + 1)))
            r_seq.append( (start, 0, (leg + 1)) )
            for ell in range(self.__len__(), len(atks)+1, -1):
                r_seq.append( (start, -1*ell, (leg + 1)))
            for (aj_r, aj_c) in reversed(atks):
                r_seq.append( (start, -abs(v_mu[aj_r]), (leg + self.leg_length(aj_r, aj_c) + 1)))
            r_seq.append( (start, 0, (leg + 1/2)) )
            r_seq = [coroot(a, b, K, self.__len__()) for (a, b, K) in r_seq]
        return(r_seq)

    def box_greedy_root_sequence(self):
        r"""
        Returns the full root sequence for the box greedy reduced work by concatenating the
        root sequences for each box in spiral order using self.root_sequence(r, c)
        """
        return(list(itertools.chain.from_iterable([self.root_sequence(r, c) for (r, c) in self.boxes()])))

    def Weyl_action(self, i):
        r"""
        returns self modified by the action of s_i, the ith simple generator of the 
        affine Weyl group of type CC_n.

        EXAMPLES:: 
        sage: D = KoornwinderBoxDiagram([0, 2, 3, -1, 1])

        sage: D.Weyl_action(0)
        [1, 2, 3, -1, 1]
        sage: D.Weyl_action(1)
        [2, 0, 3, -1, 1]
        sage: D.Weyl_action(2)
        [0, 3, 2, -1, 1]
        sage: D.Weyl_action(3)
        [0, 2, -1, 3, 1]
        sage: D.Weyl_action(4)
        [0, 2, 3, 1, -1]
        sage: D.Weyl_action(5)
        [0, 2, 3, -1, -1]

        sage: #Check that the action squares to the identity
        sage: D.Weyl_action(0).Weyl_action(0)
        [0, 2, 3, -1, 1]

        """
        if i < 0 or i > self.__len__():
            raise(IndexError("list index out of range"))
        with self.clone() as elt2:
            if i == 0:
                elt2.__setitem__(i, 1-self.__getitem__(0))
            elif i == self.__len__():
                elt2.__setitem__(i-1, -self.__getitem__(i-1))
            else:
                elt2.__setitem__(i-1, self.__getitem__(i))
                elt2.__setitem__(i, self.__getitem__(i-1))
            out = elt2
        out._compute_v_mu()
        return(elt2)
        
class AlcoveWalkTableaux(ClonableList):
    r"""
    Class for uncompressed set valued tableaux representing an alcove walk for the box greedy
    reduced word.  Initialize with shape mu and sequence of booleans for each letter in the 
    bgrw indicating whether the letter is a fold.

    EXAMPLES::

    sage: D = KoornwinderBoxDiagram([0, 2, 3, -1, 1])
    sage: bgrw = D.box_greedy_reduced_word()
    sage: z = [1, 2, 3, 4, 5]
    sage: fold_positions = [3, 6, 8, 15, 17, 22, 26, 32, 37, 38]
    sage: AW = AlcoveWalkTableaux(D, [i not in fold_positions for i in range(len(bgrw))], z)
    sage: AW.pretty_print()
                     |0|                                                                
                     | |                                                                
                     | |_____________ ___________________                               
                     |1|             |                   |                              
                     | |1, 0         |4, -, 4, 3, 2, 1, -|                              
                     | |_____________|___________________|____________________________  
                     |2|             |                   |                            | 
                     | |2, -, 0      |4, 5, 4, -, 2, 1, 0|1, 2, -, 4, 5, 4, 3, -, -, 0| 
     ________________| |_____________|___________________|____________________________| 
    |                |3|                                                                
    |5, 4, 3, 2, 1, -| |                                                                
    |________________| |_____________                                                   
                     |4|             |                                                  
                     | |4, -, 2, -, 0|                                                  
                     | |_____________|  
    sage: AW.permutation_sequence(2, 1)
    [[-2, 3, 1, 4, 5], [-2, 3, 1, 4, 5], [2, 3, 1, 4, 5]]
    sage: AW.permutation_sequence(1, 2)
    [[-4, -2, 1, 5, 3], [-4, -2, 1, 5, 3], [-4, -2, 1, 3, 5], [-4, -2, 3, 1, 5],
    [-4, 3, -2, 1, 5], [3, -4, -2, 1, 5], [3, -4, -2, 1, 5]]

    TODO: implement for alcove walks of type other that box greedy reduced word.
    """
    def __init__(self, mu, folds, z = None):
        if z is None:
            z = [i for i in range(1, len(mu)+1)]
        ClonableList.__init__(self, parent = Parent(), lst = folds)
        self.check()
        self._mu = KoornwinderBoxDiagram(mu)
        self._n = len(mu)
        self._walk_type = self._mu.box_greedy_reduced_word()
        if len(self._walk_type) != self.__len__():
            raise(ValueError("folds do not match lenght of bgrw for shape"))
        self._z = z
        self._generate_endpoints()

    def check(self):
        self._set_list([bool(tf) for tf in self._get_list()])

    def n(self):
        r"""returns the lenght of the underlying weight"""
        return(self._n)

    def shape(self):
        r"""
        get the underlying shape for self as KoornwinderBoxDiagram
        """
        return(self._mu)

    def interior_box(self, r, c):
        r"""
        Call to self._mu.interior_box to determine if (r, c) is an interior box.  
        See KoornwinderBoxDiageram.interior_box.
        """
        return(self._mu.interior_box(r, c))

    def _check_box_call(self, r, c):
        r"""
        Raises ValueError if (r, c) is not an interior box
        """
        if not self.interior_box(r, c):
            raise(ValueError("expected coordinates (r, s) to be within the shape mu"))

    def boxes(self):
        r"""
        Returns a list of boxes in underlying shape, represented as coordinate pairs 
        (r, s), and in the spiral reading order.
        """
        return(self._mu.boxes())

    def bgrw_filling(self, r, c):
        r"""
        returns the portion of the box greedy reduced word in box (r, c) of the underlying shape
        """
        return(self._mu.bgrw_filling(r, c))

    def root_sequence(self, r, c):
        r"""
        returns the portion of the root sequence in box (r, c) of the underlying shape
        """
        return(self._mu.root_sequence(r, c))

    def _generate_endpoints(self):
        r"""
        Generate (once and for all) the underlying weight for self as KoornwinderBoxDiagram
        and the permutation sequence for self as self._perm_sequence.  Accessible by 
        get_full_permutation_sequence and get_permutation_sequence.
        """
        perm_seq = [self._z]
        for i in range(self.__len__()):
            new_perm = perm_seq[-1].copy()
            b = self.__getitem__(i)
            if b:
                v = self._walk_type[i]
                if v == 0:
                    new_perm[0] = - new_perm[0]
                elif v == len(self._mu) and  b:
                    new_perm[-1] = - new_perm[-1]
                else:
                    np_i_mi = new_perm[v-1]
                    np_i = new_perm[v]
                    new_perm[v-1] = np_i
                    new_perm[v] = np_i_mi
            perm_seq.append(new_perm)
        self._perm_sequence = perm_seq

    def full_permutation_sequence(self):
        r"""
        Returns the full permutation sequences of self
        """
        return(self._perm_sequence)

    def permutation_sequence(self, r, c):
        r"""
        Returns the sequence of signed permutation associated to the box (r, s)
        """
        self._check_box_call(r, c)
        ind = self._mu.boxes().index((r, c))
        ct = 0
        first = 0
        for i in range(len(self._walk_type)):
            l = self._walk_type[i]
            if l == 0:
                if ct != ind:
                    first = i+1
                    ct += 1
                else:
                    return(self._perm_sequence[first+1:i+2])

    def folding(self, r, c):
        r"""
        Returns folding sequence, as a list of booleans, associated to the box (r, s).
        """
        self._check_box_call(r, c)
        ind = self._mu.boxes().index((r, c))
        word = self._mu.bgrw_filling(r, c)
        ct = 0
        first = 0
        for i in range(len(self._walk_type)):
            l = self._walk_type[i]
            if l == 0:
                if ct != ind:
                    first = i+1
                    ct += 1
                else:
                    return([(word[j-first], self.__getitem__(j)) for j in range(first, i+1)])

    def pretty_print(self):
        r"""
        Pretty print the argument in an ascii_art style.
        """
        lst = self._mu
        coord_dicts = {}
        for (r, c) in self._mu.boxes():
            word = ""
            for (letter, fold) in self.folding(r, c):
                if fold:
                    word += str(letter) + ", "
                else:
                    word += "-, "
            coord_dicts[(r, c)] = word[:-2]
        mx_ln = {(lambda x : x if x < 0 else x - 1)(c) : max([len(coord_dicts[k]) for k in coord_dicts if k[1] == c]) for (b, c) in coord_dicts.keys()}
        transition_intervals = []
        bar_intervals = []
        lm = 0
        rm = 0
        for i in range(len(lst)+1):
            if i == len(lst):
                transition_intervals.append(sorted([0, lst[-1]]))
                bar_intervals.append(sorted([0, lst[-1]]))
            else:
                lm = min(lst[i], lm)
                rm = max(lst[i], rm)
                if i == 0:
                    transition_intervals.append(sorted([0, lst[i]]))
                    bar_intervals.append([0, -1])
                else:
                    transition_intervals.append(sorted([lst[i-1], lst[i]]))
                    bar_intervals.append(sorted([0, lst[i-1]]))
        mx_ln[lm-1] = 1
        mx_ln[rm] = 1
        print_string = ""
        for j in range(lm, rm+1):
            ti = transition_intervals[0]
            if j == 0:
                print_string += " _"
            if ti[0] <= j and j < ti[1]:
                print_string += " " + mx_ln[j] * "_"
            else:
                print_string += " " + mx_ln[j] * " "
        print_string += "\n"
        for i in range(len(lst)):
            ti = transition_intervals[i+1]
            bi = bar_intervals[i+1]
            for j in range(lm, rm+1):
                if j == 0:
                    print_string += f"|{i}"
                if bi[0] <=j and j < bi[1]:
                    print_string += "|" + mx_ln[j] * " "
                elif j == bi[1] and j < ti[1]:
                    print_string += "|" + mx_ln[j] * " "
                elif j == bi[1] and j == ti[1]:
                    print_string += "|" + mx_ln[j] * " "
                elif ti[0] <= j and j < ti[1]:
                    print_string += " " + mx_ln[j] * " "
                else:
                    print_string += " " + mx_ln[j] * " "
            print_string += "\n"
            for j in range(lm, rm+1):
                if j == 0:
                    print_string += f"| "
                if bi[0] <=j and j < bi[1]:
                    if j >= 0:
                        k = (i, j+1)
                    else:
                        k = (i, j)
                    if k in coord_dicts:
                        print_string += "|" + coord_dicts[k] + (mx_ln[j] - len(coord_dicts[k]))*" "
                    else:
                        print_string += "|" + mx_ln[j]*" "
                elif j == bi[1] and j < ti[1]:
                    print_string += "|" + mx_ln[j] * " "
                elif j == bi[1] and j == ti[1]:
                    print_string += "|" + mx_ln[j] * " "
                elif ti[0] <= j and j < ti[1]:
                    print_string += " " + mx_ln[j] * " "
                else:
                    print_string += " " + mx_ln[j] * " "
            print_string += "\n"
            for j in range(lm, rm+1):
                if j == 0:
                    print_string += f"|_"
                if bi[0] <=j and j < bi[1]:
                    print_string += "|" + mx_ln[j] * "_"
                elif j == bi[1] and j < ti[1]:
                    print_string += "|" + mx_ln[j] * "_"
                elif j == bi[1] and j == ti[1]:
                    print_string += "|" + mx_ln[j] * " "
                elif ti[0] <= j and j < ti[1]:
                    print_string += " " + mx_ln[j] * "_"
                else:
                    print_string += " " + mx_ln[j] * " "
            print_string += "\n"
        print(print_string)

def AlcoveWalks(shape, z = None):
        r"""
        Iterable over all alcove walks for the box greedy reduced word of shape.

        EXAMPLES:: 
        sage: count = 0
        sage: for AW in AlcoveWalks([1, -1], z = None):
        ....:     count += 1
        ....:     if count % 5 == 0:
        ....:         AW.pretty_print()
        
        [True, False, True, True]
        [[1, 2], [-1, 2], [-1, 2], [2, -1], [-2, -1]]
                 _ _  
                |0| | 
                | |0| 
         _______|_|_| 
        |       |1|   
        |-, 1, 0| |   
        |_______|_|   

        [False, True, True, False]
        [[1, 2], [1, 2], [1, -2], [-2, 1], [-2, 1]]
                 _ _  
                |0| | 
                | |-| 
         _______|_|_| 
        |       |1|   
        |2, 1, -| |   
        |_______|_|   

        [False, False, False, True]
        [[1, 2], [1, 2], [1, 2], [1, 2], [-1, 2]]
                 _ _  
                |0| | 
                | |-| 
         _______|_|_| 
        |       |1|   
        |-, -, 0| |   
        |_______|_|    
        """
        mu = KoornwinderBoxDiagram(shape)
        bgrw = mu.box_greedy_reduced_word()
        if z is None:
            z = [i for i in range(1, len(mu)+1)]
        for folds in itertools.product([True, False], repeat = len(bgrw)):
            yield(AlcoveWalkTableaux(mu, folds, z = z))

class _CF_index_class(ClonableElement):
    r"""
    A small wrapper for abstract products of C and F functions, stored as a dictionary of tuples.  
    Don't actually initialize from here.
    """

    def __init__(self, parent, weights):
        ClonableElement.__init__(self, parent=parent)
        self._data = weights
        self._n = parent.n()
        self.check()
        self.set_immutable()
        
    def check(self):
        r"""
        """
        if self._data in self.parent():
            self._data = self._data.get_data()
        elif not isinstance(self._data, dict):
            raise(ValueError(f"weights passed are not of type dict (type {type(self._data)})"))
        for w in self._data:
            for function_type in self._data[w]:
                if function_type not in ["C", "F+", "F-"]:
                    raise(ValueError("expected dictionary values to specify C, F+, or F-"))
                multiplicity = self._data[w][function_type]
                if Integer(multiplicity) < 0:
                    raise(ValueError("negative powers of functions are not allowed"))
            for function_type in ["C", "F+", "F-"]:
                if function_type not in self._data[w]:
                    self._data[w][function_type] = 0
            for x in w[:-1]:
                Integer(x)
            Integer(2*w[-1])
            if len(w) != self.n() + 1:
                raise(ValueError("weights passed are of inconsistent lengtt"))

    def n(self):
        r"""
        """
        return(self._n)

    def _repr_(self):
        r"""
        """
        repr = []
        for w in self._data:
            for function_type in self._data[w]:
                exp = self._data[w][function_type]  
                if exp > 0:
                    if function_type == "C":
                       str_rep = "C_{" + str(w) + "}"
                    elif function_type == "F+":
                       str_rep = "F^{+}_{" + str(w) + "}"
                    elif function_type == "F-":
                       str_rep = "F^{-}_{" + str(w) + "}"
                    if exp > 1:
                        str_rep = "(" + str_rep + ")^{" + str(exp) + "}"
                    repr.append(str_rep)
        repr.sort()
        final_str = ""
        for str_rep in repr:
            final_str += str_rep + "*"
        final_str = final_str[:-1]
        if len(final_str) == 0:
            final_str = "1"
        return(final_str)
        
    def _hash_(self):
        r"""
        """
        return(self._repr_().__hash__())

    def __iter__(self):
        r"""
        """
        return iter(self._data)

    def get_data(self):
        r"""
        """
        return(self._data)

    def Weyl_action(self, i):
        r"""
        implements Weyl group action on tuples
        """
        if i < 0 or i > self.n():
            raise(IndexError("list index out of range"))
        new_data = {}
        for tup in self.get_data():
            new_tup = [v for v in tup]
            if i == 0:
                new_tup[-1] += new_tup[0]
                new_tup[0] *= -1
            elif i == self._n:
                new_tup[-2] *= -1
            else:
                tvi_m = new_tup[i-1]
                tvi = new_tup[i]
                new_tup[i-1] = tvi
                new_tup[i] = tvi_m
            new_tup = tuple(new_tup)
            new_data[new_tup] = self._data[tup]
        return(self.parent()._element_constructor_(new_data))

    def _add_(self, other):
        r"""
        Adds two objects by taking the union of the underlying frozen set.
        """
        if other not in self.parent():
            raise(ValueError(f"cannot add KoornwinderWeightTuples of lengths {self.n()} and {other.n()}"))
        new_elt = {w : {k : self._data[w][k] for k in self._data[w]} for w in self._data}
        for w in other._data:
            if w in new_elt:
                for function_type in other._data[w]:
                    if function_type in new_elt[w]:
                        new_elt[w][function_type] += other._data[w][function_type]
                    else:
                        new_elt[w][function_type] = other._data[w][function_type]
            else:
                new_elt[w] = {k : other._data[w][k] for k in other._data[w]}
        return(self.parent()._element_constructor_(new_elt))

class CF_indices(UniqueRepresentation, Parent):
    """
    A small parent for indices of products of C and F functions.  
    Initialize using _element_constructor_.
    """

    Element = _CF_index_class

    def __init__(self, n):
        r"""
        """
        Parent.__init__(self, category = Sets())
        self._n = n
        self._name = f"Class for products of Koornwinder C and F functions of length {n}"

    def _element_constructor_(self, weights):
        r"""
        """
        return(self.element_class(self, weights))

    def n(self):
        return(self._n)
    
    def __contains__(self, x):
        r"""
        """
        if isinstance(x, self.element_class):
            if x.n() == self._n:
                return(True)
            else:
                return(False)
   
    def _coerce_map_from_(self, S):
        try:
            if self._element_constructor_(S):
                return(True)
        except:
            pass
        return(super()._coerce_map_from_(S))


def KoornwinderPolynomialRing(n, ground_ring = None):
    r"""initializes polynomial ring for Koornwinder polynomials"""
    if ground_ring is None:
        ground_ring = QQ
    R.<sqrt_q, sqrt_t, sqrt_t0, sqrt_tn, sqrt_u0, sqrt_un> = PolynomialRing(ground_ring)
    R.latex_variable_names()[0] = '\\sqrt{q}'
    R.latex_variable_names()[1] = '\\sqrt{t}'
    R.latex_variable_names()[2] = '\\sqrt{t_{0}}'
    R.latex_variable_names()[3] = '\\sqrt{t_{n}}'
    R.latex_variable_names()[4] = '\\sqrt{u_{0}}'
    R.latex_variable_names()[5] = '\\sqrt{u_{n}}'
    CoeffRing = FractionField(R)
    KPR = LaurentPolynomialRing(CoeffRing, [f"x_{i+1}" for i in range(n)])
    return(KPR)

def _get_tu_params(KPR, wt):
    r"""
    gets t and u parameters for the coroot wt
    """
    scalars = KPR.base_ring().gens()
    if len([i for i in range(len(wt)-1) if wt[i] != 0]) == 2:  
        #print("(co)Root is in orbit 5")
        sq_t_alpha = scalars[1]
        sq_u_alpha = scalars[1]
    else:
        try:
            Integer(wt[-1])
            #print("(co)Root is in orbit 1")
            sq_t_alpha = scalars[3]
            sq_u_alpha = scalars[2]
        except: 
            #print("(co)Root is in orbit 3")
            sq_t_alpha = scalars[5]
            sq_u_alpha = scalars[4]
    #print(f"returning {(KPR(sq_t_alpha), KPR(sq_u_alpha))}")
    return(KPR(sq_t_alpha), KPR(sq_u_alpha))

def _get_Y_wt(KPR, wt, mu = None):
    r"""
    evaluates Y^wt at the character `mu' (defaulting=0)
    wt lives in Z^n \oplus Z \frac{K}{2}
    """
    n = len(wt)-1
    if mu is None:
        mu = KoornwinderBoxDiagram(n*[0])
    v_mu = mu.v_mu()
    scalars = KPR.base_ring().gens()
    sqrt_q = scalars[0]
    sqrt_t = scalars[1]
    sqrt_t0 = scalars[2]
    sqrt_tn = scalars[3]
    Y_wt = KPR.one()
    for i in range(n):
        KPR.one()
        q_front = sqrt_q^(-2*mu[i])
        t_front = sqrt_t^(-2*v_mu[i])
        inner_term = (sqrt_t0 * sqrt_tn * sqrt_t^(2*n))^(1)
        Yi_val = q_front * t_front * inner_term
        Y_wt *= Yi_val^(wt[i])
    Y_wt *= sqrt_q^(-2*wt[-1])
    #print(f"Y_{wt} evaluates to {Y_wt}")
    return(KPR(Y_wt))

def C_eval(wt, mu = None):
    r"""
    Evaluates the C function for wt at the character `mu' (defaulting=0).
    In  Macdonald's notation this would be c_wt(y^{-1}).
    Yj maps to t^(-j) and Y^{K/2} maps to q^(-1/2)
    """
    n = len(wt)-1
    if mu is None:
        mu = KoornwinderBoxDiagram(n*[0])
    KPR = KoornwinderPolynomialRing(n)
    (sq_t_alpha, sq_u_alpha) = _get_tu_params(KPR, wt)
    Y_wt = _get_Y_wt(KPR, [-x for x in wt], mu)
    num1 = (1 - sq_t_alpha * sq_u_alpha * Y_wt)
    num2 = (1 + sq_t_alpha * sq_u_alpha^(-1) * Y_wt)
    denom = (1 - Y_wt^(2))
    #print(f"{sq_t_alpha^(-1)} * ({num1}) * ({num2}) / ({denom})")
    return( sq_t_alpha^(-1) * num1 * num2 / denom)

def Fplus_eval(wt, mu = None):
    r"""
    Evaluates the F+ function for wt at the character `mu' (default=0).
    F^+_\alpha = t_\alpha^(-1/2) - C_\alpha
    """
    n = len(wt)-1
    if mu is None:
        mu = KoornwinderBoxDiagram(n*[0])
    KPR = KoornwinderPolynomialRing(n)
    (sq_t_alpha, sq_u_alpha) = _get_tu_params(KPR, wt)
    return( C_eval(wt, mu) - sq_t_alpha^(+1) )

def Fminus_eval(wt, mu = None):
    r"""
    Evaluates the F- function for wt at the character `mu' (default=0).
    F^+_\alpha = t_\alpha^(1/2) - C_\alpha
    """
    n = len(wt)-1
    if mu is None:
        mu = KoornwinderBoxDiagram(n*[0])
    KPR = KoornwinderPolynomialRing(n)
    (sq_t_alpha, sq_u_alpha) = _get_tu_params(KPR, wt)
    return( C_eval(wt, mu) - sq_t_alpha^(-1) )


class CF_function_space_element(CombinatorialFreeModule.Element):
    r"""
    wrapper for elements of CF_function_space.
    """
    def Weyl_action(self, i):
        out = self.parent().zero()
        for (m, c) in self:
            new_m = self.parent()._CF_ind(m.Weyl_action(i))
            out += c * self.parent().monomial(new_m)
        return(out)

    def expand(self):
        r"""
        evaluates each C and F to Koornwinder parameters
        """
        poylnomial_ring = self.parent().base_ring()
        out = poylnomial_ring.zero()
        for (m, c) in self:
            mon = poylnomial_ring(c)
            data = m.get_data()
            for w in data:
                for function_type in data[w]:
                    if function_type == "C":
                        mon *= C_eval(w)^data[w][function_type]
                    elif function_type == "F+":
                        mon *= Fplus_eval(w)^data[w][function_type]
                    elif function_type == "F-":
                        mon *= Fminus_eval(w)^data[w][function_type]
            out += mon
        return(out)

    def divided_difference(self, i):
        r"""
        Applies ith demazure operator del_i to self, using the `x variables' 
        from self.parent().xvars()
        """
        if i < 0 or i > self.parent().n():
            raise(ValueError(f"del_i for i={i} not implemented (must use i between 0 and {self.n()})"))
        out = self.parent().zero()
        for (m, c) in self:
            m = self.parent()._CF_ind(m)
            out += self._divided_difference_helper(i, c) * self.parent().B(m)
        return(out)

    def _divided_difference_helper(self, i, f):
        r"""
        Applies ith demazure operator del_i to polynomial (which must live in parent.base_ring()).
        """
        xv = self.parent().xvars()
        BR = self.parent().base_ring()
        try:
            f = BR(f)
        except:
            raise(NotImplementedError(f"Cannot apply divided difference to member of {type(f)}"))
        out = BR(0)
        if f == BR(1):
            return(out)
        dd = f.dict()
        if i == 0:
            q = BR(self.parent().q())
            for e in dd:
                ee = list(e)
                if ee[0] >= 0:
                    for a in range(1, ee[0]+1):
                        ep = type(e)([ee[0]-2*a] + ee[1:])
                        out += dd[e]*q^(a-1)*BR.monomial(ep)
                else:
                    for a in range(ee[0], 0):
                        ep = type(e)([ee[0]-2*a] + ee[1:])
                        out += dd[e]*q^(a-1)*BR.monomial(ep)
        elif i == self.parent().n():
            for e in dd:
                ee = list(e)
                if ee[-1] >= 0:
                    for a in range(1, ee[-1]+1):
                        ep = type(e)(ee[:-1] + [ee[-1]-2*a])
                        out += dd[e] * -1*BR.monomial(ep)
                else:
                    for a in range(ee[-1], 0):
                        ep = type(e)(ee[:-1] + [ee[-1]-2*a])
                        out += dd[e] * -1*BR.monomial(ep)
        else:
            for e in dd:
                ee = list(e)
                if e[i-1] >= e[i]:
                    for a in range(1, e[i-1]-e[i]+1):
                        #print(a)
                        ep = type(e)(ee[:i-1] + [ee[i-1]-a, ee[i]+a-1] + ee[i+1:])
                        out += dd[e] * BR.monomial(ep)
                elif e[i-1] < e[i]:
                    for a in range(1, e[i-1]-e[i]+1):
                        ep = type(e)(ee[:i-1] + [ee[i-1]-a, ee[i]+a-1] + ee[i+1:])
                        diff = dd[e] * -1 * BR.monomial(ep)
        return(self.parent().base_ring()(out))

    def tau(self, i, mu):
        r"""
        Applies tau_i operator in Koornwinder creation formula.
        """
        mu = KoornwinderBoxDiagram(mu)
        if i < 0 or i > self.parent().n():
            raise(ValueError(f"del_i for i={i} not implemented (must use i between 0 and {self.n()})"))
        if i == 0:
            return(self._tau_0_helper(mu))
        elif i == self.parent().n():
            return(self._tau_n_helper(mu))
        else:
            out = self._tau_i_helper(mu, i)
            #print(f"actually returning {out}")
            return(out)

    def _tau_0_helper(self, mu):
        r"""
        Applies tau_0 operator in Koornwinder creation formula
        """
        coeff_1_exponent = self.parent()._CF_ind({coroot(-1, 0, 1/2, self.parent().n()) : {"F-" : 1}})
        for i in mu.box_greedy_reduced_word():
            coeff_1_exponent = coeff_1_exponent.Weyl_action(i)
        coeff_1 = self.parent().B(coeff_1_exponent)
        n = self.parent().n()
        xv = self.parent().xvars()
        q = self.parent().q()
        t = self.parent().t()
        t0 = self.parent().t0()
        tn = self.parent().tn()
        u0 = self.parent().u0()
        ev_mu_Y1 = _get_Y_wt(self.parent().KPR(), [1]+n*[0], mu)
        coeff_2 = t0^(-1/2)*ev_mu_Y1*xv[0]
        coeff_3 = -t0^(-1/2)* ev_mu_Y1 * (xv[0]^3 + (u0^(-1/2)-u0^(1/2))*q^(1/2)*t0^(1/2)*xv[0]^2 - t0*q*xv[0])
        return(coeff_1 * self + coeff_2 * self + coeff_3 * self.divided_difference(0))

    def _tau_n_helper(self, mu):
        r"""
        Applies tau_0 operator in Koornwinder creation formula
        """
        coeff_1_exponent = self.parent()._CF_ind({coroot(self.parent().n(), 0, 0, self.parent().n()) : {"C" : 1}})
        for i in mu.box_greedy_reduced_word():
            coeff_1_exponent = coeff_1_exponent.Weyl_action(i)
        coeff_1 = self.parent().B(coeff_1_exponent)
        xv = self.parent().xvars()
        tn = self.parent().t0()
        un = self.parent().u0()
        coeff_2 = -tn^(-1/2)*(1 - tn^(1/2)*un^(1/2)*xv[-1])*(1 - tn^(1/2)*un^(-1/2)*xv[-1])
        return(coeff_1 * self + coeff_2 * self.divided_difference(self.parent().n()))

    def _tau_i_helper(self, mu, i):
        r"""
        Applies tau_i operator, 0<i<n, in Koornwinder creation formula
        """
        coeff_1_exponent = self.parent()._CF_ind({coroot(i, -i-1, 0, self.parent().n()) : {"C" : 1}})
        for j in mu.box_greedy_reduced_word():
            coeff_1_exponent = coeff_1_exponent.Weyl_action(j)
        coeff_1 = self.parent().B(coeff_1_exponent)
        xv = self.parent().xvars()
        t = self.parent().t()
        coeff_2 = +(t^(-1/2)*xv[i] - t^(1/2)*xv[i-1])
        return(coeff_1 * self + coeff_2 * self.divided_difference(i))


class CF_function_space(CombinatorialFreeModule):
    r"""
    A ring for Koornwinder polynomial coefficients, stored with unevaluated C and F functions.
    Use element method .expand() to express in the usual Koornwinder parameters by evaluating 
    at the character 0.
    """

    Element = CF_function_space_element

    def __init__(self, n, R = None):
        self._name = "Abstract C Function Ring"
        self._repr_option_bracket = True
        self._CF_ind = CF_indices(n)
        self._n = n
        if R == None:
            R = QQ
        self._KPR = KoornwinderPolynomialRing(n, R)
        self._scalars = self._KPR.base_ring() # the ring with Koornwinder parameters but no x's
        CombinatorialFreeModule.__init__(self, self._KPR, self._CF_ind, category=AlgebrasWithBasis(self._KPR), prefix = '', bracket=False)

    def B(self, data):
        r"""returns monomial indexed by data"""
        exp = self._CF_ind(data)
        return(self.monomial(exp))

    def product_on_basis(self, x, y):
        x = self._CF_ind._element_constructor_(x)
        y = self._CF_ind._element_constructor_(y)
        return self.monomial(x + y)
    
    def one_basis(self):
        return self._CF_ind({})
    
    def one(self):
        return self.monomial(self.one_basis()) 

    def n(self):
        """returns global_n"""
        return(self._n)

    def scalars(self):
        """returns the ring of scalars in which the Koornwinder parameters live"""
        return(self._scalars)

    def KPR(self):
        r"""
        Returns the Koornwinder polynomial ring containing the q, t, u parameters and the x-variables.
        """
        return(self._KPR)

    def q(self):
        """returns the parameter q from yq_R"""
        return(self._scalars.gens()[0]^2)

    def t(self):
        """returns the parameter t from CoefficientRing."""
        return(self._scalars.gens()[1]^2)

    def t0(self):
        """returns the parameter t0 from CoefficientRing."""
        return(self._scalars.gens()[2]^2)

    def tn(self):
        """returns the parameter tn from CoefficientRing."""
        return(self._scalars.gens()[3]^2)

    def u0(self):
        """returns the parameter u0 from CoefficientRing."""
        return(self._scalars.gens()[4]^2)

    def un(self):
        """returns the parameter un from CoefficientRing."""
        return(self._scalars.gens()[5]^2)

    def xvars(self):
        """returns set of x varriables"""
        return(self.base_ring().gens())

    def change_ring(self, R):
        r"""
        For this class the base ring should be considered the base ring of
        the underlying Koornwinder polynomial ring.  This method should returns 
        the base change (with this notion of base) of ``self`` to `R`.
        """
        if R is self.scalars().base_ring():
            return(self)
        raise NotImplementedError('the method change_ring() has not yet been implemented')

def AlcoveWalkBoxWeight(AWT, r, c):
    r"""
    Computes the weight of the box (r, c) in the alcove walk tableaux AWT
    """
    n = AWT.n()
    x_exponent = n * [0]
    Indices = CF_indices(n)
    CF_Ring = CF_function_space(n)
    CF_exponent = Indices._element_constructor_({})
    folding = AWT.folding(r, c)
    perm_seq = AWT.permutation_sequence(r, c)
    root_seq = AWT.root_sequence(r, c)
    if folding[-1][1]:
        if len(perm_seq) > 1:
            last_step = perm_seq[-2]
        else:
            last_step = [-1*perm_seq[-1][0]] + perm_seq[-1][1:]
        if last_step[0] > 0:
            x_exponent[last_step[0] - 1] += 1
        elif last_step[0] < 0:
            x_exponent[-1*last_step[0] - 1] += -1
    for i in range(len(folding)):
        fold = folding[i][1]
        step = folding[i][0]
        if not fold:
            z_bi = perm_seq[i]
            beta_i = root_seq[i]
            if step == 0:
                if z_bi[0] > 0:
                    key = "F-"
                elif z_bi[0] < 0:
                    key = "F+"
            elif step == n:
                if z_bi[n-1] > 0:
                   key = "F-"
                elif z_bi[n-1] < 0:
                    key = "F+"
            else:
                z_bi_i = z_bi[step-1]
                z_bi_ip = z_bi[step]
                if 0 < z_bi_i and z_bi_i < z_bi_ip: 
                    key = "F+"
                elif 0 < z_bi_i and z_bi_i > z_bi_ip:
                    key = "F-"
                elif z_bi_i < 0 and 0 < z_bi_ip:
                    key = "F-"
                elif z_bi_i < 0 and z_bi_i > z_bi_ip:
                    key = "F-"
                elif z_bi_ip < 0 and z_bi_i < z_bi_ip:
                    key = "F+"
            CF_exponent += Indices._element_constructor_({root_seq[i] : {key : 1}})
    return(CF_Ring.monomial(CF_exponent) * prod([CF_Ring.xvars()[i] ** x_exponent[i] for i in range(n)]))


def Weyl_len_s(zed):    
    r"""
    Computes the number of long roots ($\epsilon_i - \epsilon_j$ or $\epsilon_i + \epsilon_j$)
    negated by the signed permutation zed.
    """
    n = len(zed)
    zed_rho = [numpy.sign(zi)*(n+1-abs(zi)) for zi in zed]
    count = 0
    for i in range(len(zed_rho)):
        for j in range(i+1, len(zed_rho)):
            if zed_rho[i] + zed_rho[j] < 0:
                count += 1
            if zed_rho[i] - zed_rho[j] < 0:
                count += 1
    return(count)

def Weyl_len_d(zed):
    r"""
    Computes the number of short roots ($\epsilon_i$) negated by the signed permutation zed.
    """
    count = 0
    for i in range(len(zed)):
        if zed[i] < 0:
            count += 1
    return(count)

def AlcoveWalkWeight(AWT):
    r"""
    Given a AlcoveWalkTableaux, computes its weight as a monomial in the 
    corresponding CF_function_space
    """
    n = AWT.n()
    CF_Ring = CF_function_space(n)
    t = CF_Ring.t()
    tn = CF_Ring.tn()
    out = 1
    bx = AWT.boxes()
    if len(bx) == 0:
        return(out)
    for (r, c) in bx:
        out *= AlcoveWalkBoxWeight(AWT, r, c)
    zr = AWT.permutation_sequence(r, c)[-1]
    ls_zr = Weyl_len_s(zr)
    ld_zr = Weyl_len_d(zr)
    out *= t ** ((ls_zr )/2) * tn ** ((ld_zr)/2)
    return(out)

def KoornwinderAlcoveWalkTableau(mu, z = None):
    r"""
    Generates the relative Koornwinder polynomial E_mu^z using the tableau formulation of  
    the alcove walk formula. Returns output in terms of F functions in CF_function_space.
    """
    n = len(mu)
    mu = KoornwinderBoxDiagram(mu)
    boxes = mu.boxes()
    CF_Ring = CF_function_space(n)
    t = CF_Ring.t()
    tn = CF_Ring.tn()
    term = CF_Ring(0)
    for AW in AlcoveWalks(mu, z):
        term += AlcoveWalkWeight(AW)
    if z is None:
        z = [i+1 for i in range(n)]
    v = mu.v_mu()
    signed_perm_index = (lambda x : v.index(x) if x in v else v.index(-x))
    zv = [numpy.sign(v[signed_perm_index(i)]) * z[signed_perm_index(i)] for i in range(1, len(z)+1)]
    ls_zv = Weyl_len_s(zv)
    ld_zv = Weyl_len_d(zv)
    return(t ** (-ls_zv/2) * tn ** (-ld_zv/2) * term)

def greedy_step(mu):
    r"""
    Returns the box greedy reduced word for the KoornwinderBoxDiagram mu
    """
    n = len(mu)
    if all([mui == 0 for mui in mu]):
        return( None )
    elif all([mui <= 0 for mui in mu]):
        largest_neg = max([i for i in range(1, n+1) if mu[i-1] < 0])
        new_mu = mu.Weyl_action(largest_neg)
        return( (new_mu, largest_neg))
    else:
        smallest_pos = min([i for i in range(1, n+1) if mu[i-1] > 0])
        new_mu = mu.Weyl_action(smallest_pos-1)
        return( (new_mu, smallest_pos-1) )


def KoornwinderCreation(mu):
    r"""
    Generates the electronic Koornwinder polynomial E_mu using the creation formula.
    """
    n = len(mu)
    mu = KoornwinderBoxDiagram(mu)
    CF_Ring = CF_function_space(n)
    t = CF_Ring.t()
    tn = CF_Ring.tn()
    v = mu.v_mu()
    ls_v = Weyl_len_s(v)
    ld_v = Weyl_len_d(v)
    return(t ** (-ls_v/2) * tn ** (-ld_v/2) * _KoornwinderCreation_helper(mu))

def _KoornwinderCreation_helper(mu):
    r"""
    Uses recursive calls to the the creation formula functions, including
    - 
    - 
    to compute the un-normaled Koornwinder polynomial \hat{E}_mu.
    """
    mu = KoornwinderBoxDiagram(mu)
    CF_Ring = CF_function_space(len(mu))
    if all([v == 0 for v in mu]):
        return(CF_Ring.B({}))
    else:
        (nu, i) = greedy_step(mu)
        prev = _KoornwinderCreation_helper(nu)
        out = prev.tau(i, nu)
        return(out)

    










