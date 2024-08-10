# Pair Generator Script

This Python script generates all possible pairs from a list of entries provided by the user.

## Usage

To use the script, run the following command:

```bash
python pairwise_list.py --entries <entry1> <entry2> ... <entryN>

Replace <entry1>, <entry2>, ..., <entryN> with the entries you want to generate pairs from. For example:

python pairwise_list.py --entries apple banana cherry

This will output:

Entries: ['apple', 'banana', 'cherry']

All possible pairs: [('apple', 'banana'), ('apple', 'cherry'), ('banana', 'cherry')]

Arguments

--entries: A list of entries to pair. Provide one or more entries separated by spaces.
Example

To generate pairs from the list ["a", "b", "c"], you would run:

python pairwise_list.py --entries a b c

License

This project is licensed under the MIT License - see the LICENSE file for details.

Author

Dr. Brandon D. Jordan
