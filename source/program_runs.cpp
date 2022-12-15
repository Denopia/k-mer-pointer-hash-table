#include "program_runs.hpp"


/*
    Normal hashmap mode.
*/

int run_mode_0(int argc, char const* argv[])
{
    std::cout << "Mode 0 run starts\n";

    // --- Initialize needed arguments ---

    // k-mer length
    int k = 0;
    // hash table slots
    int slots = 100000;
    // max number of uniq entries in the hash table
    int uniq_entries_limit = 0;
    // run queries after hash map is built
    bool query_enabled = false;
    // include reversed reads
    bool reverse_reads_enabled = false;
    // first full k-mer of the read is checked before anything else
    bool first_kmer_check_enabled = false;
    // reads file path
    std::string reads_path = "No file";
    // Help trigger
    bool user_wants_help = false;
    // Verbosity flag
    bool verbose = false;
    // Only canonical flag
    bool only_canonical = false;

    int min_abundance = 1;
    bool write_output = false;
    std::string output_path = "ragavan";
    std::string output_directory = ".";

    uint64_t genome_length_estimate = 0;
    uint64_t read_coverage_estimate = 0;

    int argi = 1;

    // --- Parse arguments ---
    while (argi < argc)
    {
        std::string as(argv[argi]);
        if (as.compare("-m") == 0){
            // skip mode argument
            argi += 2;
        }
        else if (as.compare("-k") == 0){
            k = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-s") == 0)
        {
            slots = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-u") == 0)
        {
            uniq_entries_limit = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-q") == 0)
        {
            query_enabled = true;
            argi += 1;
        }
        else if (as.compare("-r") == 0)
        {
            reverse_reads_enabled = true;
            argi += 1;
        }
        else if (as.compare("-c") == 0)
        {
            only_canonical = true;
            argi += 1;
        }
        else if (as.compare("-f") == 0)
        {
            first_kmer_check_enabled = true;
            argi += 1;
        }
        else if (as.compare("-p") == 0)
        {
            reads_path = std::string(argv[argi+1]);
            argi += 2;
        }
        else if (as.compare("-h") == 0)
        {
            user_wants_help = true;
            argi += 1;
        }
        else if (as.compare("-v") == 0)
        {
            verbose = true;
            argi += 1;
        }
        else if (as.compare("-a") == 0){
            min_abundance = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-o") == 0)
        {
            output_path = std::string(argv[argi+1]);
            write_output = true;
            argi += 2;
        }
        else if (as.compare("-d") == 0)
        {
            output_directory = std::string(argv[argi+1]);
            argi += 2;
        }
    }
    // At the moment 8 different arguments accepted + help

    if (user_wants_help)
    {   std::cout << "This is a k-mer has table program. Here are the accepted arguments:\n\n";
        std::cout << "\t-m [integer] determines the program mode\n";
        std::cout << "\t-k [integer] determines the k-mer length\n";
        std::cout << "\t-s [integer] determines the minimum hash table size (actual used size is the next prime)\n";
        std::cout << "\t-u [integer] determines a limit for the number of items that can be stored in the hash table (default: maximum size)\n";
        std::cout << "\t-q flag enables queries after hash table is built\n";
        std::cout << "\t-r flag tells the program to also use the reversed reads to build the hash table\n";
        std::cout << "\t-f flag tells the program to do a first k-mer check after buffer is empty (should reduce the number of stored entries)\n";
        std::cout << "\t-p [string] indicates the path to the read file (at the moment fasta format is required)\n";
        std::cout << "\t-v flag tells the program to print what is going on during the program\n";
        std::cout << "\t-c flag tells the program to only consider canonical k-mers\n";
        std::cout << "\t-h flag shows this help message\n\n";

        return 0;
    }

    // --- Next, verify some arguments ---    
     if (reads_path == "No File")
    {
        std::cout << "** ERROR ** I need a path to the reads file (must be fasta)...\n";
        return 1;
    }
    if (k < 1)
    {
        std::cout << "** ERROR ** The value of k must be greater than 0...\n";
        return 1;
    }
    // cout vai oliko joku cerr?

    // --- Make the number of hash table slots prime ---
    uint32_t prime_slots = mathfunctions::next_prime(slots);
    
    // --- Define new variable for hash table size ---
    int hash_table_slots = prime_slots;

    float load_factor = 0.8;

    int resizing_count = 0;



    // --- Set item limit to max if not specified by user ---
    //if (item_limit == 0)
    //    item_limit = prime_slots;        
    
    // --- Build file reader ---
    FastaReader* file_reader = new FastaReader(reads_path, reverse_reads_enabled);

    // --- Build k-mer factory ---
    KMerFactoryCanonical2BC* kmer_factory = new KMerFactoryCanonical2BC(k);

    KMerFactoryCanonical2BC* resizing_kmer_factory = new KMerFactoryCanonical2BC(k);
    
    // --- Build rolling hasher ----
    RollingHasher1* hasher = new RollingHasher1(hash_table_slots, k);
    //AdderHasher1* hasher = new AdderHasher1(hash_table_slots);

    // --- Build hash table ---
    BasicHashTable1* hash_table = new BasicHashTable1(hash_table_slots, k);

    /*

        TODO:
            1) Check fasta reader is ok - done
            2) Check k-mer factory is ok -
            3) Check rolling hasher is ok -
            4) Check hash table is ok -
            5) hpp -> cpp constructor ??

    */

    
    uint64_t current_kmer_slot;
    uint64_t previous_kmer_slot;

    uint64_t current_kmer_hash = 0;

    //bool previous_kmer_exists;

    int read_position;
    //int first_check_read_position;
    bool first_kmer_unread;

    int reads_processed = 0;

    // --- Starting reading the reads and adding found k-mers to the hash table ---
    while (file_reader->read_is_loaded())
    {
        

        // if read is too short it is skipped
        if (file_reader->get_current_read_length() < k)
        {
            file_reader->roll_to_next_read();
            continue;
        }

        read_position = 0;

        // Start reading characters from the read until it is fully read
        while (read_position < file_reader->get_current_read_length())
        {

            // Push character into the k-mer factory buffer
            kmer_factory->push_new_character(file_reader->get_current_read_character_at(read_position));
            // If it was not a legit character, reset hahser too
            
            if (kmer_factory->get_number_of_stored_characters() == 0){
                hasher->reset();
            } else {
                hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
            }

            //std::cout << "File reader character is " << file_reader->get_current_read_character_at(read_position) << "\n";
            //std::cout << "Newest character is "  << kmer_factory->get_newest_character() << "\n";
            //std::cout << "Oldest character is "  << kmer_factory->get_pushed_off_character() << "\n";
            //std::cout << "k-mer hash is " << hasher->get_current_hash() << "\n\n";

            // If k-mer factory has k characters stored, update the hash table
            if (kmer_factory->get_number_of_stored_characters() == k)
            {
                if (only_canonical)
                {
                    if (kmer_factory->forward_kmer_is_canonical())
                        hash_table->find_and_increment_andifnot_insert(kmer_factory, hasher->get_current_hash());
                }
                else
                {
                    hash_table->find_and_increment_andifnot_insert(kmer_factory, hasher->get_current_hash());
                }
            }

            if (hash_table->get_number_of_inserted_items() > load_factor*hash_table_slots)
            {
                resizing_count += 1;

                std::cout << " * Resizing happening at read " << reads_processed << " position " << read_position-1 << "\n";
                //std::cout << "Last character pushed was: " << file_reader->get_current_read_character_at(read_position-1) << "\n";
                std::cout << " * Resizing hash table\n";

                std::cout << " * Estimating new size\n";

                // If estimates are not given, just double
                
                if (genome_length_estimate == 0 || read_coverage_estimate == 0)
                {
                    prime_slots = mathfunctions::next_prime(2*hash_table_slots);
                }
                else
                {
                    if (hash_table->get_number_of_inserted_items() < genome_length_estimate)
                    {
                        prime_slots = mathfunctions::next_prime(2*hash_table_slots);
                    }
                    else
                    {
                        prime_slots = mathfunctions::next_prime(2*hash_table_slots);
                    }

                }
                
                hash_table_slots = prime_slots;
                
                std::cout << " * Creating new objects for the new hash table\n";

                // FIND THE SIZE OF THE NEW HASH TABLE
                // --- Build a new rolling hasher ----
                //delete hasher;
                RollingHasher1 * new_hasher = new RollingHasher1(hash_table_slots, k);
                
                delete hasher;
                hasher = new_hasher;
                new_hasher = NULL;

                // --- Build a new hash table ---
                BasicHashTable1* hash_table_new = new BasicHashTable1(hash_table_slots, k);

                // COPY CONTENT FROM THE OLD TABLE INTO THE NEW ONE
                std::cout << " * Start copying k-mers in to the new hash table\n";
                resizing_kmer_factory->reset();
                hash_table_new->copy_content_from_another_hash_table(resizing_kmer_factory, hasher, hash_table);
                resizing_kmer_factory->reset();

                delete hash_table;
                hash_table = hash_table_new;
                hash_table_new = NULL;

                hasher->reset();
                hasher->load_full_factory_canonical(kmer_factory);
                std::cout << "===============================================================================\n";

            }
            read_position+=1;
        }

        file_reader->roll_to_next_read();
        kmer_factory->reset();
        hasher->reset();

        reads_processed += 1;

        if (verbose && (reads_processed % 100 == 0))
        {
            std::cout << "This was read number " << reads_processed << "\n";
            std::cout << hash_table->get_number_of_inserted_complete_kmers() << "/" << hash_table_slots << " k-mers stored currently\n";
            std::cout << "====================================================================\n";
        }
    }

    if (verbose)
    {
        std::cout << "This was read number " << reads_processed << "\n";
        std::cout << hash_table->get_number_of_inserted_complete_kmers() << "/" << hash_table_slots << " k-mers stored currently\n";
        std::cout << "====================================================================\n";
    }

    if (write_output)
    {
        std::chrono::high_resolution_clock::time_point start_time, end_time;

        start_time = std::chrono::high_resolution_clock::now();

        std::cout << "Writing k-mers to output file\n";
        hash_table->write_kmers_on_disk(min_abundance, output_path);

        end_time = std::chrono::high_resolution_clock::now();

        uint64_t writing_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end_time-start_time).count();

        std::cout << "### Writing time: " << writing_microseconds << " microseconds\n";
    }

    // --- Run queries if they are enabled ---
    if (query_enabled)
    {
        std::string query_kmer;
        if (verbose)
            std::cout << "Query k-mers (input q to quit):\n";
        while (true)
        {
            current_kmer_slot = 0;

            std::cin >> query_kmer;
            if (!query_kmer.compare("q") || !query_kmer.compare("Q"))
            {
                if (verbose)
                    std::cout << "Query end now.\n";
                break;
            }    
            else if (query_kmer.length()!=uint64_t(k))
                if (verbose)
                    std::cout << "k-mer must be " << k << " long...\n";
                else
                    std::cout << "-1\n";
            else if (!twobitstringfunctions::is_clean_string(query_kmer))
                if (verbose)
                    std::cout << "Your query k-mer is not valid...\n";
                else
                    std::cout << "-1\n";
            else 
            {
                kmer_factory->reset();
                hasher->reset();

                for (char cc : query_kmer)
                {
                    kmer_factory->push_new_character(cc);
                    hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                }
                current_kmer_hash = hasher->get_current_hash();
                current_kmer_slot = hash_table->find(kmer_factory, current_kmer_hash);

                if (current_kmer_slot == hash_table_slots){
                    std::cout << "0\n";
                } else {
                    std::cout << hash_table->get_kmer_count_in_slot(current_kmer_slot) << "\n";
                }
            }
        }
    }

    // --- Free allocated memory ---
    delete file_reader;
    delete kmer_factory;
    delete hasher;
    // Hash table was not deleted in previous version ?? 
    delete hash_table; 

    // --- Program run is now finished ---
    std::cout << "Total number of solid (count at least 2): " << hash_table->get_number_of_solid_kmers() << "\n";
    std::cout << "Mode 0 run ends\n";

    return 0;

}


/*
    Single character with pointer to previous k-mer hashmap mode.
*/
int run_mode_1(int argc, char const* argv[])
{
    std::cout << "Mode 1 run starts\n";

    // --- Initialize needed arguments ---

    // k-mer length
    int k = 0;
    // hash table slots
    int slots = 100000;
    // max number of uniq entries in the hash table
    int uniq_entries_limit = 0;
    // run queries after hash map is built
    bool query_enabled = false;
    // include reversed reads
    bool reverse_reads_enabled = false;
    // first full k-mer of the read is checked before anything else
    bool first_kmer_check_enabled = false;
    // reads file path
    std::string reads_path = "No file";
    // Help trigger
    bool user_wants_help = false;
    // Verbosity flag
    bool verbose = false;
    // Only canonical flag
    bool only_canonical = false;

    int min_abundance = 1;
    bool write_output = false;
    std::string output_directory = ".";

    int argi = 1;

    // --- Parse arguments ---
    while (argi < argc)
    {
        std::string as(argv[argi]);
        if (as.compare("-m") == 0){
            // skip mode argument
            argi += 2;
        }
        else if (as.compare("-k") == 0){
            k = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-s") == 0)
        {
            slots = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-u") == 0)
        {
            uniq_entries_limit = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-q") == 0)
        {
            query_enabled = true;
            argi += 1;
        }
        else if (as.compare("-r") == 0)
        {
            reverse_reads_enabled = true;
            argi += 1;
        }
        else if (as.compare("-c") == 0)
        {
            only_canonical = true;
            argi += 1;
        }
        else if (as.compare("-f") == 0)
        {
            first_kmer_check_enabled = true;
            argi += 1;
        }
        else if (as.compare("-p") == 0)
        {
            reads_path = std::string(argv[argi+1]);
            argi += 2;
        }
        else if (as.compare("-h") == 0)
        {
            user_wants_help = true;
            argi += 1;
        }
        else if (as.compare("-v") == 0)
        {
            verbose = true;
            argi += 1;
        }
        else if (as.compare("-a") == 0){
            min_abundance = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-o") == 0)
        {
            write_output = true;
            argi += 1;
        }
        else if (as.compare("-d") == 0)
        {
            output_directory = std::string(argv[argi+1]);
            argi += 2;
        }
        
    }
    // At the moment 8 different arguments accepted + help

    if (user_wants_help)
    {   std::cout << "This is a k-mer has table program. Here are the accepted arguments:\n\n";
        std::cout << "\t-m [integer] determines the program mode\n";
        std::cout << "\t-k [integer] determines the k-mer length\n";
        std::cout << "\t-s [integer] determines the minimum hash table size (actual used size is the next prime)\n";
        std::cout << "\t-u [integer] determines a limit for the number of items that can be stored in the hash table (default: maximum size)\n";
        std::cout << "\t-a [integer] minimum k-mer abundance (for writing output, default: 1)\n";
        std::cout << "\t-q flag enables queries after hash table is built\n";
        std::cout << "\t-r flag tells the program to also use the reversed reads to build the hash table\n";
        std::cout << "\t-f flag tells the program to do a first k-mer check after buffer is empty (should reduce the number of stored entries)\n";
        std::cout << "\t-p [string] indicates the path to the read file (at the moment fasta format is required)\n";
        std::cout << "\t-v flag tells the program to print what is going on during the program\n";
        std::cout << "\t-c flag tells the program to only consider canonical k-mers\n";
        std::cout << "\t-o flag tell the program to write found k-mers into output file\n";
        std::cout << "\t-d [string] shows the directory where output is written\n";
        std::cout << "\t-h flag shows this help message\n\n";

        return 0;
    }

    // --- Next, verify some arguments ---    
     if (reads_path == "No File")
    {
        std::cout << "** ERROR ** I need a path to the reads file (must be fasta)...\n";
        return 1;
    }
    if (k < 1)
    {
        std::cout << "** ERROR ** The value of k must be greater than 0...\n";
        return 1;
    }
    // cout vai oliko joku cerr?

    // --- Make the number of hash table slots prime ---
    uint32_t prime_slots = mathfunctions::next_prime(slots);
    
    // --- Define new variable for hash table size ---
    int hash_table_slots = prime_slots;

    // --- Set item limit to max if not specified by user ---
    //if (item_limit == 0)
    //    item_limit = prime_slots;        
    
    // --- Build file reader ---
    FastaReader* file_reader = new FastaReader(reads_path, reverse_reads_enabled);

    // --- Build k-mer factory ---
    KMerFactory2BC* kmer_factory = new KMerFactory2BC(k);
    
    // --- Build rolling hasher ----
    RollingHasher1* hasher = new RollingHasher1(hash_table_slots, k);

    // --- Build hash table ---
    PointerHashTable1* hash_table = new PointerHashTable1(hash_table_slots, k);

    /*

        TODO:
            1) Check fasta reader is ok - done
            2) Check k-mer factory is ok -
            3) Check rolling hasher is ok -
            4) Check hash table is ok -
            5) hpp -> cpp constructor ??

    */

    uint64_t current_kmer_slot = 0;
    uint64_t previous_kmer_slot = hash_table_slots;

    uint64_t current_kmer_hash = 0;

    //bool previous_kmer_exists = false;

    int read_position = 0;
    int first_check_read_position = 0;
    bool first_kmer_unread = true;

    int reads_processed = 0;

    // --- Starting reading the reads and adding found k-mers to the hash table ---
    while (file_reader->read_is_loaded())
    {
        
        // if read is too short it is skipped
        if (file_reader->get_current_read_length() < k)
        {
            file_reader->roll_to_next_read();
            reads_processed += 1;
            continue;
        }

        read_position = 0;
        previous_kmer_slot = hash_table_slots;

        // Start reading characters from the read until it is fully read
        while (read_position < file_reader->get_current_read_length())
        {
            // if k-mer factory is empty, we can perform first k-mer check
            // also, the first k-mer check needs to be enabled
            if ((kmer_factory->get_number_of_stored_characters() == 0) && (first_kmer_check_enabled))
            {
                // First, set the position where we start the first check
                first_check_read_position = read_position;
                // Then we add characters to the k-mer factory until it is full OR read is fully read
                while ((kmer_factory->get_number_of_stored_characters() < k) && (first_check_read_position < file_reader->get_current_read_length()))
                {
                    // Push character into the k-mer factory buffer
                    kmer_factory->push_new_character(file_reader->get_current_read_character_at(first_check_read_position));
                    // If it was not a legit character, reset hahser too
                    if (kmer_factory->get_number_of_stored_characters() == 0){
                        hasher->reset();
                    } else {
                        hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                    }
                    first_check_read_position+=1;
                }
                // Once k-mer factory is full or read is fully read, we look if the k-mer is already in the hash table or move to the next read
                //std::cout << "ENOUGH CHARACTERS STORED\n";
                if (kmer_factory->get_number_of_stored_characters() == k)
                {
                    current_kmer_slot = hash_table->find_and_increment(kmer_factory, hasher->get_current_hash(), hash_table_slots);
                    // If k-mer does not exist, we need to go back and add incomplete k-mers
                    if (current_kmer_slot == hash_table_slots)
                    {
                        //previous_kmer_exists = false;
                        kmer_factory->reset();
                        hasher->reset();
                    }
                    // Otherwise, we can continue from here
                    else
                    {
                        //previous_kmer_exists = true;
                        previous_kmer_slot = current_kmer_slot;
                        read_position = first_check_read_position;
                    }
                }
                // If read has no more characters, break out
                else
                {
                    read_position = first_check_read_position;
                }
            }
            // Check if read positions is still legal, kinda annoying that this needs to be done...
            if (read_position >= file_reader->get_current_read_length())
            {
                break;
            }

            // Push character into the k-mer factory buffer
            kmer_factory->push_new_character(file_reader->get_current_read_character_at(read_position));
            // If it was not a legit character, reset hahser too
            if (kmer_factory->get_number_of_stored_characters() == 0){
                hasher->reset();
            } else {
                hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
            }

            //std::cout << "File reader character is " << file_reader->get_current_read_character_at(read_position) << "\n";
            //std::cout << "Newest character is "  << kmer_factory->get_newest_character() << "\n";
            //std::cout << "Oldest character is "  << kmer_factory->get_pushed_off_character() << "\n";
            //std::cout << "k-mer hash is " << hasher->get_current_hash() << "\n\n";

            // If k-mer factory has k characters stored, update the hash table
            if (kmer_factory->get_number_of_stored_characters() > 0)
            {
                current_kmer_slot = hash_table->find_and_increment_andifnot_insert(kmer_factory, hasher->get_current_hash(), previous_kmer_slot);
            }
            read_position+=1;
            previous_kmer_slot = current_kmer_slot;
        }

        file_reader->roll_to_next_read();
        kmer_factory->reset();
        hasher->reset();

        reads_processed += 1;

        if (verbose)
        {
            std::cout << "This was read number " << reads_processed << "\n";
            std::cout << hash_table->get_number_of_inserted_complete_kmers() << "/" << hash_table_slots << " k-mers stored currently\n";
            std::cout << hash_table->get_number_of_inserted_items() << "/" << hash_table_slots << " items  stored currently\n";
            std::cout << "====================================================================\n";
        }
        

    }

    // --- Run queries if they are enabled ---
    if (query_enabled)
    {
        std::string query_kmer;
        if (verbose)
            std::cout << "Query k-mers (input q to quit):\n";
        while (true)
        {
            current_kmer_slot = 0;

            std::cin >> query_kmer;
            if (!query_kmer.compare("q") || !query_kmer.compare("Q"))
            {
                if (verbose)
                    std::cout << "Query end now.\n";
                break;
            }    
            else if (query_kmer.length()!=uint64_t(k))
                if (verbose)
                    std::cout << "k-mer must be " << k << " long...\n";
                else
                    std::cout << "-1\n";
            else if (!twobitstringfunctions::is_clean_string(query_kmer))
                if (verbose)
                    std::cout << "Your query k-mer is not valid...\n";
                else
                    std::cout << "-1\n";
            else 
            {
                kmer_factory->reset();
                hasher->reset();

                for (char cc : query_kmer)
                {
                    kmer_factory->push_new_character(cc);
                    hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                }
                current_kmer_hash = hasher->get_current_hash();
                current_kmer_slot = hash_table->find(kmer_factory, current_kmer_hash);

                if (current_kmer_slot == hash_table_slots){
                    std::cout << "0\n";
                } else {
                    std::cout << hash_table->get_kmer_count_in_slot(current_kmer_slot) << "\n";
                }
            }
        }
    }

    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    int countM = 0;
    for (uint32_t l = 0; l < hash_table_slots; l++){
        if (hash_table->kmer_in_slot_is_complete(l)){
            if (hash_table->get_kmer_count_in_slot(l) == 1){
                count1 += 1;
            } else if (hash_table->get_kmer_count_in_slot(l) == 2){
                count2 += 1;
            } else if (hash_table->get_kmer_count_in_slot(l) == 3){
                count3 += 1;
            } else if (hash_table->get_kmer_count_in_slot(l) > 3){
                countM += 1;
            }
        }
    }

    std::cout << "k-mer abundance = 1 : " << count1 << "\n";
    std::cout << "k-mer abundance = 2 : " << count2 << "\n";
    std::cout << "k-mer abundance = 3 : " << count3 << "\n";
    std::cout << "k-mer abundance > 3 : " << countM << "\n";


    // --- Free allocated memory ---
    delete file_reader;
    delete kmer_factory;
    delete hasher;
    // Hash table was not deleted in previous version ?? 
    delete hash_table; 

    // --- Program run is now finished ---
    std::cout << "Total number of solid (count at least 2): " << hash_table->get_number_of_solid_kmers() << "\n";
    std::cout << "Mode 1 run ends\n";

    return 0;

}

/*
    Single character with pointer to previous k-mer hashmap mode.
    Temp storage to avoid inserting incomplete k-mers
*/
int run_mode_2(int argc, char const* argv[])
{
    std::cout << "Mode 2 run starts\n";

    // --- Initialize needed arguments ---

    // k-mer length
    int k = 0;
    // hash table slots
    int slots = 0;
    // max number of uniq entries in the hash table
    int uniq_entries_limit = 0;
    // run queries after hash map is built
    bool query_enabled = false;
    // include reversed reads
    bool reverse_reads_enabled = false;
    // first full k-mer of the read is checked before anything else
    bool first_kmer_check_enabled = false;
    // reads file path
    std::string reads_path = "No file";
    // Help trigger
    bool user_wants_help = false;
    // Verbosity flag
    bool verbose = false;
    // Only canonical flag
    bool only_canonical = false;

    int min_abundance = 1;
    bool write_output = false;
    std::string output_path = "ragavan";
    std::string output_directory = ".";

    uint64_t genome_length_estimate = 0;
    uint64_t read_coverage_estimate = 0;

    float estimated_error_rate = 0.01;

    int resizing_count = 0;

    int argi = 1;

    // --- Parse arguments ---
    while (argi < argc)
    {
        std::string as(argv[argi]);
        if (as.compare("-m") == 0){
            // skip mode argument
            argi += 2;
        }
        else if (as.compare("-k") == 0){
            k = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-s") == 0)
        {
            slots = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-u") == 0)
        {
            uniq_entries_limit = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-q") == 0)
        {
            query_enabled = true;
            argi += 1;
        }
        else if (as.compare("-r") == 0)
        {
            reverse_reads_enabled = true;
            argi += 1;
        }
        else if (as.compare("-c") == 0)
        {
            only_canonical = true;
            argi += 1;
        }
        else if (as.compare("-f") == 0)
        {
            first_kmer_check_enabled = true;
            argi += 1;
        }
        else if (as.compare("-p") == 0)
        {
            reads_path = std::string(argv[argi+1]);
            argi += 2;
        }
        else if (as.compare("-h") == 0)
        {
            user_wants_help = true;
            argi += 1;
        }
        else if (as.compare("-v") == 0)
        {
            verbose = true;
            argi += 1;
        }
        else if (as.compare("-a") == 0){
            min_abundance = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-g") == 0){
            genome_length_estimate = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-x") == 0){
            read_coverage_estimate = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-o") == 0)
        {
            output_path = std::string(argv[argi+1]);
            write_output = true;
            argi += 2;
        }
        else if (as.compare("-d") == 0)
        {
            output_directory = std::string(argv[argi+1]);
            argi += 2;
        }
    }
    // At the moment 8 different arguments accepted + help

    if (user_wants_help)
    {   std::cout << "This is a k-mer has table program. Here are the accepted arguments:\n\n";
        std::cout << "\t-m [integer] determines the program mode\n";
        std::cout << "\t-k [integer] determines the k-mer length\n";
        std::cout << "\t-a [integer] minimum k-mer abundance (for writing output, default: 1)\n";
        std::cout << "\t-s [integer] determines the minimum hash table size (actual used size is the next prime)\n";
        std::cout << "\t-u [integer] determines a limit for the number of items that can be stored in the hash table (default: maximum size)\n";
        std::cout << "\t-q flag enables queries after hash table is built\n";
        std::cout << "\t-r flag tells the program to also use the reversed reads to build the hash table\n";
        std::cout << "\t-f flag tells the program to do a first k-mer check after buffer is empty (should reduce the number of stored entries)\n";
        std::cout << "\t-p [string] indicates the path to the read file (at the moment fasta format is required)\n";
        std::cout << "\t-v flag tells the program to print what is going on during the program\n";
        std::cout << "\t-c flag tells the program to only consider canonical k-mers\n";
        std::cout << "\t-o flag tell the program to write found k-mers into output file\n";
        std::cout << "\t-d [string] shows the directory where output is written\n";
        std::cout << "\t-g [integer] estimated length of the genome\n";
        std::cout << "\t-x [integer] estimated read coverage\n";
        std::cout << "\t-h flag shows this help message\n\n";

        return 0;
    }

    // --- Next, verify some arguments ---    
     if (reads_path == "No File")
    {
        std::cout << "** ERROR ** I need a path to the reads file (must be fasta)...\n";
        return 1;
    }
    if (k < 1)
    {
        std::cout << "** ERROR ** The value of k must be greater than 0...\n";
        return 1;
    }
    if (slots == 0)
    {
        std::cout << "Initial number of hash table slots not given, trying to estimate using genome length and coverage instead\n";
        if ((genome_length_estimate == 0) || (read_coverage_estimate == 0))
        {
            std::cout << "** ERROR ** If hash table slots value is not given, genome length estimate and read coverage must be given\n";
            exit(1);
        }
        else
        {
            float incorrect_kmers_portion = 1.0 - std::pow((1.0-estimated_error_rate), k);
            slots = std::ceil(incorrect_kmers_portion * genome_length_estimate * read_coverage_estimate) + genome_length_estimate;
        }
    }
    // cout vai oliko joku cerr?

    // --- Make the number of hash table slots prime ---
    uint32_t prime_slots = mathfunctions::next_prime(slots);
    
    // --- Define new variable for hash table size ---
    int hash_table_slots = prime_slots;

    float load_factor = 0.8;

    // --- Set item limit to max if not specified by user ---
    //if (item_limit == 0)
    //    item_limit = prime_slots;        
    
    // --- Build file reader ---
    FastaReader* file_reader = new FastaReader(reads_path, reverse_reads_enabled);

    // --- Build k-mer factory ---
    KMerFactory2BC* kmer_factory = new KMerFactory2BC(k);

     // --- Build second k-mer factory for resizing---
    KMerFactory2BC* resizing_kmer_factory = new KMerFactory2BC(k);
    
    // --- Build rolling hasher ----
    RollingHasher1* hasher = new RollingHasher1(hash_table_slots, k);

    // --- Build hash table ---
    PointerHashTable2* hash_table = new PointerHashTable2(hash_table_slots, k, kmer_factory->number_of_blocks);

    /*

        TODO:
            1) Check fasta reader is ok - done
            2) Check k-mer factory is ok -
            3) Check rolling hasher is ok -
            4) Check hash table is ok -
            5) hpp -> cpp constructor ??

    */

    uint64_t current_kmer_slot = 0;
    uint64_t previous_kmer_slot = hash_table_slots;

    uint64_t current_kmer_hash = 0;

    bool previous_kmer_exists = false;

    int read_position = 0;
    int first_check_read_position = 0;
    bool first_kmer_unread = true;

    int reads_processed = 0;

    uint64_t counter_10k = 0;
    uint64_t counter_100k = 0;
    uint64_t counter_1m = 0;
    uint64_t counter_10m = 0;

    uint64_t added_kmers_at_last_10k_mark = 0;
    uint64_t added_kmers_at_last_100k_mark = 0;
    uint64_t added_kmers_at_last_1m_mark = 0;
    uint64_t added_kmers_at_last_10m_mark = 0;

    uint64_t last_10k_interval_new_kmers = 0;
    uint64_t last_100k_interval_new_kmers = 0;
    uint64_t last_1m_interval_new_kmers = 0;
    uint64_t last_10m_interval_new_kmers = 0;


    // --- Starting reading the reads and adding found k-mers to the hash table ---
    while (file_reader->read_is_loaded())
    {
        
        // if read is too short it is skipped
        if (file_reader->get_current_read_length() < k)
        {
            file_reader->roll_to_next_read();
            reads_processed += 1;
            continue;
        }

        previous_kmer_exists = false;
        previous_kmer_slot = hash_table_slots;
        read_position = 0;

        // Start reading characters from the read until it is fully read
        while (read_position < file_reader->get_current_read_length())
        {
            // Push character into the k-mer factory buffer
            kmer_factory->push_new_character(file_reader->get_current_read_character_at(read_position));
            // If it was not a legit character, reset hahser too
            if (kmer_factory->get_number_of_stored_characters() == 0){
                hasher->reset();
                previous_kmer_exists = false;
            } else {
                hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
            }

            //std::cout << "File reader character is " << file_reader->get_current_read_character_at(read_position) << "\n";
            //std::cout << "Newest character is "  << kmer_factory->get_newest_character() << "\n";
            //std::cout << "Oldest character is "  << kmer_factory->get_pushed_off_character() << "\n";
            //std::cout << "k-mer hash is " << hasher->get_current_hash() << "\n\n";

            // If k-mer factory has k characters stored, update the hash table
            if (kmer_factory->get_number_of_stored_characters() == k)
            {
                current_kmer_slot = hash_table->find_and_increment_andifnot_insert(kmer_factory, hasher->get_current_hash(), previous_kmer_exists, previous_kmer_slot, false); 
                previous_kmer_exists = true;
                counter_10k += 1;
                counter_100k += 1;
                counter_1m += 1;
                counter_10m += 1;

                // Update counters to track how many new k-mers have ben added in the past
                if (counter_10k == 10000)
                {
                    last_10k_interval_new_kmers = hash_table->get_number_of_inserted_items() - added_kmers_at_last_10k_mark;
                    added_kmers_at_last_10k_mark = hash_table->get_number_of_inserted_items();
                    counter_10k = 0;
                    if (counter_100k == 100000)
                    {
                        last_100k_interval_new_kmers = hash_table->get_number_of_inserted_items() - added_kmers_at_last_100k_mark;
                        added_kmers_at_last_100k_mark = hash_table->get_number_of_inserted_items();
                        counter_100k = 0;
                        if (counter_1m == 1000000)
                        {
                            last_1m_interval_new_kmers = hash_table->get_number_of_inserted_items() - added_kmers_at_last_1m_mark;
                            added_kmers_at_last_1m_mark = hash_table->get_number_of_inserted_items();
                            counter_1m = 0;
                            if (counter_10m == 10000000)
                            {
                                last_10m_interval_new_kmers = hash_table->get_number_of_inserted_items() - added_kmers_at_last_10m_mark;
                                added_kmers_at_last_10m_mark = hash_table->get_number_of_inserted_items();
                                counter_10m = 0;
                            }
                        }
                    }
                }
            }
            read_position+=1;
            previous_kmer_slot = current_kmer_slot;

            // If hash table is too full, resize
            if (hash_table->get_number_of_inserted_items() > load_factor*hash_table_slots)
            {
                resizing_count += 1;

                std::cout << " * Resizing happening at read " << reads_processed << " position " << read_position-1 << "\n";
                //std::cout << "Last character pushed was: " << file_reader->get_current_read_character_at(read_position-1) << "\n";
                std::cout << " * Resizing hash table\n";

                std::cout << " * Estimating new size\n";

                // If estimates are not given, just double
                
                if (genome_length_estimate == 0 || read_coverage_estimate == 0)
                {
                    prime_slots = mathfunctions::next_prime(2*hash_table_slots);
                }
                else
                {
                    if (hash_table->get_number_of_inserted_items() < genome_length_estimate)
                    {
                        prime_slots = mathfunctions::next_prime(2*hash_table_slots);
                    }
                    else
                    {
                        prime_slots = mathfunctions::next_prime(2*hash_table_slots);
                    }

                }
                
                //prime_slots = mathfunctions::next_prime(2*hash_table_slots);
                hash_table_slots = prime_slots;
                
                // FIRST MIGRATE ALL TEMP K-MERS IN TO THE MAIN HASH TABLE
                //std::cout << "*** CLEANING HASH TABLE TEMP ARRAY ***\n";
                //hash_table->clean_up_temp_array(kmer_factory, hasher);

                std::cout << " * Creating new objects for the new hash table\n";

                // FIND THE SIZE OF THE NEW HASH TABLE
                // --- Build a new rolling hasher ----
                //delete hasher;
                RollingHasher1 * new_hasher = new RollingHasher1(hash_table_slots, k);
                
                delete hasher;
                hasher = new_hasher;
                new_hasher = NULL;

                // --- Build a new hash table ---
                PointerHashTable2* hash_table_new = new PointerHashTable2(hash_table_slots, k, resizing_kmer_factory->number_of_blocks);

                // COPY CONTENT FROM THE OLD TABLE INTO THE NEW ONE
                std::cout << " * Start copying k-mers in to the new hash table\n";
                resizing_kmer_factory->reset();
                hash_table_new->copy_content_from_another_hash_table(resizing_kmer_factory, hasher, hash_table);
                resizing_kmer_factory->reset();

                delete hash_table;
                hash_table = hash_table_new;
                hash_table_new = NULL;

                hasher->reset();
                hasher->load_full_factory(kmer_factory);

                previous_kmer_slot = hash_table->find(kmer_factory, hasher->get_current_hash());
                if (previous_kmer_slot >= hash_table_slots)
                {
                    std::cout << "WEIRD ERROR!!!!!!!\n";
                }

                std::cout << "===============================================================================\n";

            }

        }

        file_reader->roll_to_next_read();
        kmer_factory->reset();
        hasher->reset();

        reads_processed += 1;

        if (verbose && (reads_processed % 100 == 0))
        {
            std::cout << "Processed read number " << reads_processed << "\n";
            std::cout << hash_table->get_number_of_inserted_complete_kmers() << "/" << hash_table_slots << " k-mers stored currently in MAIN\n";
            std::cout << hash_table->get_number_of_inserted_items() << "/" << hash_table_slots << " items  stored currently in MAIN\n";
            std::cout << "--------------------------------------------------------\n";
            std::cout << hash_table->get_number_of_max_temp_slots_in_use() << "/" << hash_table->get_number_of_max_temp_slots() << " k-mer slots touched in TEMP\n";
            std::cout << hash_table->get_number_of_temp_slots_in_use() << "/" << hash_table->get_number_of_max_temp_slots() << " k-mers stored currently in TEMP\n";
            std::cout << "--------------------------------------------------------\n";
            std::cout << "This many new k-mers were added in the latest 10k k-mer interval:  " << last_10k_interval_new_kmers / 10000.0 << "\n";
            std::cout << "This many new k-mers were added in the latest 100k k-mer interval: " << last_100k_interval_new_kmers / 100000.0 << "\n";
            std::cout << "This many new k-mers were added in the latest 1m k-mer interval:   " << last_1m_interval_new_kmers / 1000000.0 << "\n";
            std::cout << "This many new k-mers were added in the latest 10m k-mer interval:  " << last_10m_interval_new_kmers / 10000000.0 << "\n";
            std::cout << "===============================================================================\n";
        }
    }

    std::cout << "Migrating k-mers from TEMP to MAIN array\n";

    hash_table->clean_up_temp_array(kmer_factory, hasher);
    //hash_table->clean_up_temp_array_OLD(kmer_factory, hasher);


    std::cout << hash_table->get_number_of_inserted_complete_kmers() << "/" << hash_table_slots << " k-mers stored currently in MAIN\n";
    std::cout << hash_table->get_number_of_inserted_items() << "/" << hash_table_slots << " items  stored currently in MAIN\n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << hash_table->get_number_of_max_temp_slots_in_use() << "/" << hash_table->get_number_of_max_temp_slots() << " k-mer slots touched in TEMP\n";
    std::cout << hash_table->get_number_of_temp_slots_in_use() << "/" << hash_table->get_number_of_max_temp_slots() << " k-mers stored currently in TEMP\n";
    std::cout << "=========================================================\n";

    if (write_output)
    {
        std::chrono::high_resolution_clock::time_point start_time, end_time;

        start_time = std::chrono::high_resolution_clock::now();

        std::cout << "Writing k-mers to output file\n";
        hash_table->write_kmers_on_disk(kmer_factory, min_abundance, only_canonical, output_path);

        end_time = std::chrono::high_resolution_clock::now();

        uint64_t writing_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(end_time-start_time).count();

        std::cout << "### Writing time: " << writing_microseconds << " microseconds\n";
    }


    // --- Run queries if they are enabled ---
    if (query_enabled)
    {
        std::string query_kmer;
        if (verbose)
            std::cout << "Query k-mers (input q to quit):\n";
        while (true)
        {
            current_kmer_slot = 0;

            std::cin >> query_kmer;
            if (!query_kmer.compare("q") || !query_kmer.compare("Q"))
            {
                if (verbose)
                    std::cout << "Query end now.\n";
                break;
            }    
            else if (query_kmer.length()!=uint64_t(k))
                if (verbose)
                    std::cout << "k-mer must be " << k << " long...\n";
                else
                    std::cout << "-1\n";
            else if (!twobitstringfunctions::is_clean_string(query_kmer))
                if (verbose)
                    std::cout << "Your query k-mer is not valid...\n";
                else
                    std::cout << "-1\n";
            else 
            {
                kmer_factory->reset();
                hasher->reset();

                for (char cc : query_kmer)
                {
                    kmer_factory->push_new_character(cc);
                    hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                }
                current_kmer_hash = hasher->get_current_hash();
                current_kmer_slot = hash_table->find(kmer_factory, current_kmer_hash);

                if (current_kmer_slot == hash_table_slots){
                    std::cout << "0\n";
                } else {
                    std::cout << hash_table->get_kmer_count_in_slot(current_kmer_slot) << "\n";
                }
            }
        }
    }

    std::cout << "@@@ Hash table had to be resized " << resizing_count << " times @@@\n";

    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    int countM = 0;
    for (uint32_t l = 0; l < hash_table_slots; l++){
        if (hash_table->kmer_in_slot_is_complete(l)){
            if (hash_table->get_kmer_count_in_slot(l) == 1){
                count1 += 1;
            } else if (hash_table->get_kmer_count_in_slot(l) == 2){
                count2 += 1;
            } else if (hash_table->get_kmer_count_in_slot(l) == 3){
                count3 += 1;
            } else if (hash_table->get_kmer_count_in_slot(l) > 3){
                countM += 1;
            }
        }
    }


    std::cout << "k-mer abundance = 1 : " << count1 << "\n";
    std::cout << "k-mer abundance = 2 : " << count2 << "\n";
    std::cout << "k-mer abundance = 3 : " << count3 << "\n";
    std::cout << "k-mer abundance > 3 : " << countM << "\n";

    std::cout << "- CUMULATIVE COUNTS -\n";
    std::cout << "k-mer abundance >= 1 : " << count1+count2+count3+countM << "\n";
    std::cout << "k-mer abundance >= 2 : " << count2+count3+countM << "\n";
    std::cout << "k-mer abundance >= 3 : " << count3+countM << "\n";
    std::cout << "k-mer abundance >  3 : " << countM << "\n";

    // --- Free allocated memory ---
    //delete file_reader;
    //delete kmer_factory;
    //delete hasher;
    // Hash table was not deleted in previous version ?? 
    //delete hash_table; 

    // --- Program run is now finished ---
    std::cout << "Total number of solid (count at least 2): " << hash_table->get_number_of_solid_kmers() << "\n";
    std::cout << "Mode 2 run ends\n";

    return 0;

}
