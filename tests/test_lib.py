import unittest
import tempfile
from source.lib import get_stable_configs

class TestLib(unittest.TestCase):
    def setUp(self):
        test_tbn_string = "6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)"

        self.tbn_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
        self.tbn_filename = self.tbn_file.name
        self.tbn_file.write(test_tbn_string)
        self.tbn_file.close()

    def test_get_stable_configs(self):
        configurations = list(get_stable_configs(self.tbn_filename))
        self.assertEqual(3, len(configurations))  # (6 + 3 - 5 - 4), (6 - 1 - 5), (6 - 2 - 4)
        self.assertEqual(2, configurations[0].number_of_polymers())
        self.assertEqual(2, configurations[1].number_of_polymers())
        self.assertEqual(2, configurations[2].number_of_polymers())
