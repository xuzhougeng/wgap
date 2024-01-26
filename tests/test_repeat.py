import unittest

from wgap.scripts.repeat import repeat_masker_loader
from wgap.scripts.repeat import Repeat

class RepeatAdapterTestCase(unittest.TestCase):
    
    def test_repeat_masker_loader(self):
        repeat_file = "tests/data/test_repeat_masker.out"
        repeats = repeat_masker_loader(repeat_file)
        repeat = repeats[0]
        # test class type
        self.assertIsInstance(repeat, Repeat)
        # test attributes
        self.assertEqual(repeat.chrom, "ptg000006l_1")

if __name__ == '__main__':
    unittest.main()