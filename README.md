# k-mer-pointer-hash-table

Uses sdsl-lite, original can be found here https://github.com/simongog/sdsl-lite

## Install

- mkdir build
- cd build 
- cmake -S ../source/ -B .
- cmake --build .

If everything goes like it should, you will find KMerHashtable executable in the build directory. (I would not be surprised if this did not work though)

## Use

./KMerHashtable [Options] -k [K] -s [S] -m [M] -p [P]

### Required 
- -k [K] = k-mer length
- -s [S] = initial hash table slots
- -m [M] = hash table mode, 0 = basic hash table, (1 = don't use this), 2 = pointer hash table
- -p [P] = path to reads file (only fasta accepted at the moment)

### Optional 
- -c = consider only canonical k-mers, used by mode 0 (should not affect mode 2)
- -v = print messages during program run (some messages may still be given even if not enabled... might be fixed sometime in the future)
- -r = read reverse complements are taken into account (this option should be used, otherwise the results won't be correct)
- -f = fast check (this option should also be used since it makes the program faster, I should make this be enabled automatcally)
- some other options I'm not listing here because they might not work correctly or don't do anything at the moment

