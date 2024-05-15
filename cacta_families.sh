#!/bin/bash
#
# Divides CACTA element sequences into families.
# Families here are distinguished based on clustering performed by VSEARCH [1].
# Thereafter clusters are filtered based on the minimum required size (default = 2).
#
# [1] - Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584

usage() {
  echo
  echo "Usage: $0 -i <input file> [-m <number>] [-e] [-g <GFF3 file>] [-c] [-h]" 1>&2
  echo
  exit 1
}

print_help() {
  echo
  echo "cacta_families.sh clusters CACTA transposons into families."
  echo "Clustering is performed by VSEARCH and clusters are filtered "
  echo "based on the minimum required size, 2 by default."
  echo "Tool generates FASTA file with elements, GFF3 annotation of "
  echo "filtered elements (GFF3 from detect_cacta.py must be provided), or "
  echo "FASTA file containing consensus sequence for each family."
  echo
  echo "Usage: $0 -i <input file> [-m <number>] [-e] [-g <GFF3 file>] [-c] [-h]"
  echo
  echo "Available options"
  echo "-i <input file>  Input FASTA file with CACTA TE sequences."
  echo "-m <number>      Minimum cluster members. Clusters with number of"
  echo "                   members lower than <number> are discarded."
  echo "                   Default value is 2."
  echo "-c               FASTA file \"<input file>_consensi.fasta\" containing"
  echo "                   consensus sequence for each family is generated."
  echo "-e               FASTA file \"<input file>_filtered.fasta\" containing"
  echo "                   filtered elements is generated."
  echo "-g <GFF3>        GFF3 annotation \"<input file>_filtered.gff3\" of "
  echo "                   filtered elements is generated. GFF3 file generated"
  echo "                   by detect_cacta.py must be provided."
  echo "-h               Prints this help message."
  echo
}

check_vsearch_installed() {
  if ! command -v vsearch >/dev/null; then
    echo
    echo "VSEARCH not found. Make sure to include it in the PATH"
    echo
    exit 1
  fi
}

parse_args() {
  while getopts ":i:m:g:che" opt; do
    case ${opt} in
    i)
      input=$OPTARG
      input_name=$(basename "$input")
      input_basename="${input_name%.*}"
      ;;
    m)
      min_size=$OPTARG
      ;;
    c)
      consensus=true
      ;;
    e)
      elements=true
      ;;
    g)
      annotation=true
      gff=$OPTARG
      gff_name=$(basename "$gff")
      gff_basename="${gff_name%.*}"
      ;;
    h)
      print_help
      exit 0
      ;;
    :)
      echo
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      usage
      ;;
    ?)
      echo
      echo "Invalid option: $OPTARG" 1>&2
      usage
      ;;
    esac
  done

  shift $((OPTIND - 1))

  if [ -z "$input" ]; then
    echo
    echo "Input file -i argument is required"
    usage
  fi

  if [ -z "$min_size" ]; then
    min_size=2
  fi
}

prepare_tmp_dir() {
  temp_dir=$(mktemp -d)

  if [[ ! "$temp_dir" || ! -d "$temp_dir" ]]; then
    # mktemp attempts to create a temporary directory in $temp_dir if specified, /tmp otherwise
    echo
    echo "Unable to create temporary directory to process results of VSEARCH"
    exit 1
  fi

  mkdir "$temp_dir"/clusters
}

function remove_tmp() {
  rm -rf "$temp_dir"
}

create_clusters() {
  local vsearch_opts consensus_path

  echo "VSEARCH clustering started"
  echo

  vsearch_opts=(
    --cluster_fast "$input"
    --clusterout_id
    --clusterout_sort
    --clusters "$temp_dir/clusters/FAM"
    --strand both
    --id 0.8
    --iddef 1
    --notrunclabels
  )

  if [ "$consensus" = true ]; then
    consensus_path="$temp_dir"/consensi.fasta
    vsearch_opts+=(--consout "$consensus_path")
  fi

  if ! vsearch "${vsearch_opts[@]}"; then
    exit 1
  fi
}

filter_clusters() {
  echo
  echo "Removing clusters that have less than $min_size members"

  for f in "$temp_dir/clusters/"*; do
    [ "$(grep -c ">" "$f")" -lt "$min_size" ] && rm "$f"
  done

  if [[ $(find "$temp_dir/clusters/" -type f | wc -l) -le 0 ]]; then
    echo
    echo "No family having at least $min_size members. No output generated"
    echo
    exit 0
  fi

  for f in "$temp_dir/clusters/"*; do
    FAMILY_NAME=$(basename "$f")
    sed -i "/^>/ s/$/_$FAMILY_NAME/" "$f"
  done

  echo "Filtering of CACTA elements finished"
  echo
}

generate_elements_file() {
  echo
  echo "Generating filtered elements FASTA file"

  cat "$temp_dir/clusters/"* >"$input_basename"_filtered.fasta
  echo >>"$input_basename"_filtered.fasta

  echo "Filtered CACTA elements stored in ""$input_basename""_filtered.fasta"
  echo
}

generate_gff_file() {
  echo
  echo "Generating filtered annotation GFF3 file"

  if ! cp "$gff" "$temp_dir/"; then
    echo "Unable to read $gff GFF3 file"
    return 1
  fi

  grep ">" "$temp_dir/clusters/"* --no-filename | cut -c2- | while IFS= read -r line; do
    replace="${line/_FAM*/}"
    sed -i "s/$replace/$line/" "$temp_dir/""$gff_name"
  done

  grep ">" "$temp_dir/clusters/"* --no-filename | cut -c2- | grep -Ff - "$temp_dir/""$gff_name" >"$gff_basename"_filtered.gff3
  echo >>"$gff_basename"_filtered.gff3

  echo "Filtered CACTA annotation stored in ""$gff_basename""_filtered.gff3"
  echo
}

generate_consensi_file() {
  echo
  echo "Generating FASTA consensus sequences file for CACTA families"

  while IFS= read -r line; do
    if [[ $line == ">"* ]]; then
      seqs=$(echo "$line" | sed -n 's/.*seqs=\([0-9]*\);.*/\1/p')

      if [[ $seqs -lt $min_size ]]; then
        break
      fi

    fi
    echo "$line"
  done <"$temp_dir"/consensi.fasta >"$input_basename"_consensi.fasta
  echo >>"$input_basename"_consensi.fasta

  echo "Family consensus sequences stored in ""$input_basename""_consensi.fasta"
  echo
}

main() {
  trap remove_tmp EXIT

  check_vsearch_installed
  parse_args "$@"
  prepare_tmp_dir
  create_clusters
  filter_clusters

  if [ "$elements" = true ]; then
    generate_elements_file
  fi

  if [ "$annotation" = true ]; then
    generate_gff_file
  fi

  if [ "$consensus" = true ]; then
    generate_consensi_file
  fi

  remove_tmp

}

main "$@"
