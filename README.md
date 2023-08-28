[![MIT License][license-shield]][license-url]


<!-- TABLE OF CONTENTS -->

<h3>Table of Contents</h3>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
      </ul>
    </li>
    <li><a href="#references-download">References Download</li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>




<!-- ABOUT THE PROJECT -->
## About The Project

The following snakefile and supporting scripts are designed to process FASTQ files from CUT&Tag. Briefly they will perform the following:

* Initial quality check with FastQC and MultiQC.
* Alignment with Bowtie2. Alignment parameters are based on the recommendation by Henikoff et.al.
* File conversion, sorting, and indexing with samtools.
* Removal of low quality reads (MAPQ<30).
* Removal of reads mapped to blacklist region.
* MACS2 peak calling.
* 500k bin count matrix formation.
* Reads normalization with TMM (normalization method from edgeR).
* Coverage files generation (bigwig files)

Steps involving quality filtering and normalization is done to correct for read depth and signal-to-noise ratio bias.

!! IMPORTANT !!  
Reads duplicate removal with Picard is not done in this workflow as CUT&Tag duplicates are likely to be real reads.

<!-- GETTING STARTED -->
## Getting Started

Clone the repository to your local folder with:

```
git clone https://github.com/ZainulArifin1/CUTandTag-Primary-Analysis.git

or with SSH key

git clone git@github.com:ZainulArifin1/CUTandTag-Primary-Analysis.git
```

## Prerequisites

To ensure a consistent and isolated environment for your bioinformatics project, you will need to install the required packages and libraries using either Conda or Mamba package managers. While both options are viable, I recommend using Mamba due to its superior speed and reliability.

#### Using Mamba (Recommended)

Mamba is a faster and more efficient alternative to Conda. Follow these steps to create and activate your environment using Mamba:

1. Open your Linux terminal or WSL for Windows user.

2. Navigate to the directory containing the project's environment.yml file.

3. Run the following command to create the environment:

```
mamba env create -f environment.yml
mamba activate bioinformatics
```

#### Using Conda (Alternative)

If you prefer to use Conda, you can achieve the same environment setup using these steps:

1. Open your Linux terminal or WSL for Windows user.

2. Navigate to the directory containing the project's environment.yml file.

3. Run the following command to create the environment:

```
conda env create -f environment.yml
conda activate bioinformatics
```

#### Important packages used in this workflow:

* bowtie2=2.5.1
* fastqc=0.12.1
* multiqc=1.15
* deeptools=3.5.1
* bedtools=2.31.0
* subread=2.0.6
* macs2=2.2.9.1
* samtools=1.17
* bioconductor-edger=3.40.0
* r-tidyverse=1.3.2
* snakemake=6.0.5
* snakemake-minimal=6.0.5

Grab a tea or whatever you want while waiting because this going to take a while.

<p align="center">
<img src="https://github.com/ZainulArifin1/CUTandTag-Primary-Analysis/blob/master/img/kermit-the-frog-sip.gif">
</p>

<!-- DOWNLOAD REFERENCE -->
## References Download

Please note that the GitHub repository does not include the indexed hg38 reference file required for your bioinformatics analysis. You will need to download and prepare this reference file separately. Follow the steps below to obtain and set up the reference files:

```
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/GRCh38_no_alt.zip
```

Do not forget to unzip the folder. Put the indexed references inside the folder "GRCh38_noalt_as".

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

If you have your FASTQ files ready, you must make sure first your folder is in the following format (**VERY IMPORTANT**).

<p align="center">
<img src="https://github.com/ZainulArifin1/CUTandTag-Primary-Analysis/blob/master/img/Folder_Order.PNG">
</p>

You can see the schematic below:

```
├── data
│   ├── blacklist
│   │   └── hg38-blacklist.v2.bed
│   ├── raw_fastq
│   │   ├── SRAEXAMPLE_R1.fastq
│   │   └── SRAEXAMPLE_R2.fastq
│   └── reference
│       ├── GRCh38_noalt_as
│       │   ├── GRCh38_noalt_as.1.bt2
│       │   ├── GRCh38_noalt_as.2.bt2
│       │   ├── GRCh38_noalt_as.3.bt2
│       │   ├── GRCh38_noalt_as.4.bt2
│       │   ├── GRCh38_noalt_as.rev.1.bt2
│       │   └── GRCh38_noalt_as.rev.2.bt2
│       ├── hg38_chrom_sizes
│       └── hg38_chrom_sizes_binned_500k.bed
├── environment.yml
├── job.s
├── LICENSE
├── README.md
├── scripts
│   ├── count_effective_genome_size.sh
│   ├── count_reads_in_bam.sh
│   ├── mergeCount.R
│   ├── remove_dup_column.R
│   ├── search_files_run64.sh
│   └── TMM.R
└── snakefile_k27me3_k4
```

Before running the snakefile, make sure of the following

1. In the example snakefile, "data/h3k27me3_h3k4/raw_fastq/" is the directory of the raw FASTQ files. Please adjust **EVERY INSTANCES** accordingly. You can do find and replace to safe time.

2. Please adjust the naming of FASTQ files. Example: "data/h3k27me3_h3k4/raw_fastq/{sra}_1.fastq". In this case the forward read is denoted as "_1" and the file extension is .fastq (can be fastq.gz). This instance is in rule bowtie2 (line 72)

3. Perform dry run (test run) with the following command

```
snakemake -np -s <snakefile_name>
```

Thats it! You have done all the hardwork and now you can just run following code and check in a few hours (or days depending on your data and resources).

```
snakemake --cores <num_of_cores> -s <snakefile_name>
```

If you are running the program in cluster, you can modify the sbatch file and run it with:

```
sbatch job.s
```

If there is an error, do not be afraid! Check a file called **errLog** and it will tell you where the error is. Should you require any help please raise an issue on GitHub or contact me through my email (muhammad.arifin@ucdconnect.ie).

<p align="center">
<img src="https://github.com/ZainulArifin1/CUTandTag-Primary-Analysis/blob/master/img/end.gif">
</p>

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ROADMAP -->
## Roadmap

- [x] Commit the workflow to GitHub
- [x] Full pipeline with minimum reproducible workflow

To customize the Snakefile for the CUT&Tag workflow, you can make adjustments based on the CUT&Tag tutorial available in this [website](https://yezhengstat.github.io/CUTTag_tutorial/).


<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Muhammad Zainul Arifin  
PhD Student, University College Dublin

* Twitter: [@SaintZainn](https://twitter.com/SaintZainn)
* Linkedin: [Muhammad Zainul Arifin](https://www.linkedin.com/in/muhammad-zainul-a-479aa1151/)
* Email: muhammad.arifin@ucdconnect.ie

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgment

I am grateful to the [Dey Lab](https://deylab.com/members/) for affording me the opportunity to contribute to the project that has led to the establishment of this GitHub repository.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/othneildrew/Best-README-Template.svg?style=for-the-badge
[contributors-url]: https://github.com/othneildrew/Best-README-Template/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/othneildrew/Best-README-Template.svg?style=for-the-badge
[forks-url]: https://github.com/othneildrew/Best-README-Template/network/members
[stars-shield]: https://img.shields.io/github/stars/othneildrew/Best-README-Template.svg?style=for-the-badge
[stars-url]: https://github.com/othneildrew/Best-README-Template/stargazers
[issues-shield]: https://img.shields.io/github/issues/othneildrew/Best-README-Template.svg?style=for-the-badge
[issues-url]: https://github.com/othneildrew/Best-README-Template/issues
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=for-the-badge
[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/othneildrew
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 