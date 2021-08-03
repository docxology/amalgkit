import numpy as np
from amalgkit.util import *
import glob


def list_duplicates(seq):
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in seq if x in seen or seen_add(x))
    # turn the set into a list (as requested)
    return list(seen_twice)


def parse_metadata(args, metadata):
    print("Checking essential entries from metadata file.")
    species = metadata.df.loc[:, 'scientific_name']
    uni_species = np.unique(species)
    if len(uni_species):
        print(len(uni_species), " species detected:")
        print(uni_species)
    else:
        raise ValueError(len(uni_species), " species detected. Please check if ", args.metadata,
                         " has a 'scientific_name' column.")

    sra_ids = metadata.df.loc[:, 'run']
    if len(sra_ids):
        print(len(np.unique(sra_ids)), " SRA runs detected:")
        print(np.unique(sra_ids))
        # check for duplicate runs
        if len(sra_ids) > len(np.unique(sra_ids)):
            raise ValueError("Duplicate SRA IDs detected, where IDs should be unique. Please check these entries: ",
                             list_duplicates(sra_ids))
    else:
        raise ValueError(len(np.unique(sra_ids)), " SRA runs detected. Please check if ", args.metadata,
                         " has a 'run' column.")
    return uni_species, sra_ids


def check_getfastq_outputs(args, sra_ids, metadata, output_dir):
    print("checking for getfastq outputs: ")
    if args.getfastq_dir:
        getfastq_path = args.getfastq
    else:
        getfastq_path = os.path.join(args.out_dir, "getfastq")
    data_available = []
    data_unavailable = []
    metadata_unavailable = []
    if os.path.exists(getfastq_path):
        print("amalgkit getfastq output folder detected. Checking presence of output files.")
        for sra_id in sra_ids:
            print("\n")
            print("Looking for ", sra_id)
            sra_path = os.path.join(getfastq_path, sra_id)
            if os.path.exists(sra_path):
                sra_stat = get_sra_stat(sra_id, metadata)
                try:
                    ext = get_newest_intermediate_file_extension(sra_stat, sra_path)
                except AssertionError:
                    print("could not find any fastq files for ", sra_id,
                          "Please make sure amalgkit getfastq ran properly")
                    data_unavailable.append(sra_id)
                    continue

                files = glob.glob(os.path.join(sra_path, sra_id + "*" + ext))
                print("found: ", files)
                data_available.append(sra_id)
                if args.updated_metadata_dir:
                    updated_metadata_path = os.path.join(args.updated_metadata_dir, "metadata_" + sra_id + ".tsv")
                else:
                    updated_metadata_path = os.path.join(args.out_dir, "metadata/updated_metadata/metadata_"
                                                         + sra_id + ".tsv")

                print("checking for updated metadata in: ", updated_metadata_path)
                if os.path.exists(updated_metadata_path):
                    print("found updated metadata!")
                else:
                    print("could not find updated metadata! Please make sure amalgkit getfastq ran properly for ",
                          sra_id)
                    metadata_unavailable.append(sra_id)
            else:
                print("Could not find getfastq output for: ", sra_id, "\n")
                print("Suggested command for rerun: getfastq -e email@adress.com --id ", sra_id, " -w ", args.out_dir, "--redo yes --gcp yes --aws yes --ncbi yes")
                data_unavailable.append(sra_id)

    else:
        print("Could not find getfastq output folder ", getfastq_path, ". Have you run getfastq yet?")
        data_unavailable = metadata.df['run'].tolist()

    if data_unavailable:
        print("writing SRA IDs without getfastq output to: ", os.path.join(output_dir, "SRA_IDs_without_fastq.txt"))
        file = open(os.path.join(output_dir, "SRA_IDs_without_fastq.txt"), "w")
        for sra_id in data_unavailable:
            file.write(sra_id + "\n")
        file.close()
    else:
        print("Sequences found for all SRA IDs in ", args.metadata, " !")

    if metadata_unavailable:
        print("writing SRA IDs without updated metadata output to: ", os.path.join(output_dir,
                                                                                   "SRA_IDs_without_metadata.txt"))
        file = open(os.path.join(output_dir, "SRA_IDs_without_metadata.txt"), "w")
        for sra_id in metadata_unavailable:
            file.write(sra_id + "\n")
        file.close()

    return data_available, data_unavailable


def check_quant_index(args, uni_species, output_dir):
    if args.index_dir:
        index_dir_path = args.index_dir
    else:
        index_dir_path = os.path.join(args.out_dir, "Index")
    index_unavailable = []
    index_available = []
    if os.path.exists(index_dir_path):
        for species in uni_species:
            sci_name = species.replace(" ", "_")
            sci_name = sci_name.replace(".", "")
            index_path = os.path.join(index_dir_path, sci_name + "*")
            print("\n")
            print("Looking for Index file", index_path, "for species: ", species)
            index_files = glob.glob(index_path)
            if not index_files:
                print("could not find anything in", index_path)
                sci_name = species.split(" ")
                # Deprecate subspecies or variants and look again
                # I.e. if Gorilla_gorilla_gorilla.idx was not found, we look for Gorilla_gorilla.idx instead.
                if len(sci_name) > 2:
                    sci_name = sci_name[0] + "_" + sci_name[1]
                    print("Ignoring subspecies.")
                    index_path = os.path.join(index_dir_path, sci_name + "*")
                    print("Looking for ", index_path)
                    index_files = glob.glob(index_path)
                    if index_files:
                        print("Found ", index_files, "!")
                        index_available.append(species)
                    else:
                        print("Could not find any index files for ", species)
                        index_unavailable.append(species)
                else:
                    print("Could not find any index files for ", species)
                    index_unavailable.append(species)

            else:
                print("Found ", index_files, "!")
                index_available.append(species)

            if len(index_files) > 1:
                print("Multiple possible Index files detected for ", species, ": ", index_files,
                      ". You may have to resolve ambiguity")

        if index_unavailable:
            print("writing species without Index to: ", os.path.join(output_dir, "species_without_index.txt"))
            file = open(os.path.join(output_dir, "species_without_index.txt"), "w")
            for species in index_unavailable:
                file.write(species + "\n")
            file.close()
        else:
            print("Index found for all species in ", args.metadata, " !")
    else:
        print("Could not find Index directory ", index_dir_path, " . Did you provide the correct Path?")

    return index_available, index_unavailable


def check_quant_output(args, sra_ids, output_dir):
    print("checking for quant outputs: ")
    quant_path = os.path.join(args.out_dir, "quant")
    data_available = []
    data_unavailable = []

    if os.path.exists(quant_path):
        print("amalgkit quant output folder detected. Checking presence of output files.")
        for sra_id in sra_ids:
            warned = []
            print("\n")
            print("Looking for ", sra_id)
            sra_path = os.path.join(quant_path, sra_id)
            if os.path.exists(sra_path):
                print("Found output folder ", sra_path, " for ", sra_id)
                print("Checking for output files.")
                abundance_file = os.path.join(sra_path, sra_id + "_abundance.tsv")
                run_info_file = os.path.join(sra_path, sra_id + "_run_info.json")
                abundance_h_file = os.path.join(sra_path, sra_id + "_abundance.h5")
                if os.path.exists(abundance_file) and os.path.exists(run_info_file) and os.path.exists(
                        abundance_h_file):
                    print("All quant output files present for, ", sra_id, " !")
                    data_available.append(sra_id)
                    continue
                elif not os.path.exists(abundance_file):
                    print(abundance_file, " is missing! Please check if quant ran correctly")
                    warned = True
                elif not os.path.exists(run_info_file):
                    print(run_info_file, " is missing! Please check if quant ran correctly")
                    warned = True
                elif not os.path.exists(abundance_h_file):
                    print(abundance_h_file, " is missing! Please check if quant ran correctly")
                    warned = True

                if warned:
                    data_unavailable.append(sra_id)

            else:
                print("Could not find output folder ", sra_path, " for ", sra_id)
                data_unavailable.append(sra_id)
    else:
        print("Could not find quant output folder ", quant_path, ". Have you run quant yet?")

    if data_unavailable:
        print("writing SRA IDs without quant output to: ", os.path.join(output_dir, "SRA_IDs_without_quant.txt"))
        file = open(os.path.join(output_dir, "SRA_IDs_without_quant.txt"), "w")
        for sra_id in data_unavailable:
            file.write(sra_id + "\n")
        file.close()
    else:
        print("Quant outputs found for all SRA IDs in ", args.metadata, " !")

    return data_available, data_unavailable


def sanity_main(args):
    print("reading metadata from:", args.metadata)
    metadata = load_metadata(args)

    output_dir = os.path.join(args.out_dir, 'sanity')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    uni_species, sra_ids = parse_metadata(args, metadata)

    if args.getfastq or args.all:
        check_getfastq_outputs(args, sra_ids, metadata, output_dir)

    if args.index or args.all:
        check_quant_index(args, uni_species, output_dir)

    if args.quant or args.all:
        check_quant_output(args, sra_ids, output_dir)
