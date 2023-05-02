/// Evidence that a given cell is iNKT or MAIT.   Each ExactSubclonotype has one
/// instantiation of this structure for iNKT and one for MAIT.
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct InvariantTCellAnnotation {
    #[prost(bool, required, tag = "1")]
    pub alpha_chain_gene_match: bool,
    #[prost(bool, required, tag = "2")]
    pub alpha_chain_junction_match: bool,
    #[prost(bool, required, tag = "3")]
    pub beta_chain_gene_match: bool,
    #[prost(bool, required, tag = "4")]
    pub beta_chain_junction_match: bool,
}
/// Representation of an alignment
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct Alignment {
    /// Start of the alignment in the reference
    #[prost(uint32, required, tag = "1")]
    pub ref_start: u32,
    /// Cigar string
    #[prost(string, required, tag = "2")]
    pub cigar: ::prost::alloc::string::String,
}
/// Defines a chain within an exact subclonotype.
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct ExactSubClonotypeChain {
    /// Nucleotide sequence of the chain. This will only contain ACGT alphabets
    #[prost(bytes = "vec", required, tag = "1")]
    pub nt_sequence: ::prost::alloc::vec::Vec<u8>,
    /// Amino acid sequence from the start codon at the beginning of the V-REGION.
    /// This can be inferred from the `nt_sequence` and `v_start`, but stored for
    /// convenience
    #[prost(bytes = "vec", required, tag = "2")]
    pub aa_sequence: ::prost::alloc::vec::Vec<u8>,
    /// Index of the start of the V-REGION in the `nt_sequence`.
    #[prost(uint32, required, tag = "3")]
    pub v_start: u32,
    /// Index of the end of the J-REGION in the `nt_sequence` (exclusive).
    #[prost(uint32, required, tag = "4")]
    pub j_end: u32,
    /// Index of the C-REGION of this chain in the universal reference.
    /// TODO: Should we store UVDJ regions here for convenience? The reason why it
    /// is not stored at this level is because all exact subclonotypes share the
    /// same UVDJ regions.
    #[prost(uint32, optional, tag = "5")]
    pub c_region_idx: ::core::option::Option<u32>,
    /// Index of the start of the CDR3 sequence in the `nt_sequence`. The start of
    /// the CDR3 amino acid in the `aa_sequence` is `(cdr3_start - v_start)/3`.
    #[prost(uint32, required, tag = "6")]
    pub cdr3_start: u32,
    /// Index of the end of the CDR3 sequence in the `nt_sequence` (exclusive).
    /// The end of the CDR3 amino acid in the `aa_sequence` is
    /// `(cdr3_end - v_start)/3`.
    #[prost(uint32, required, tag = "7")]
    pub cdr3_end: u32,
    /// UMI counts of contigs associated with this exact subclonotype chain. The
    /// number of elements in this vector is equal to the number of barcodes
    /// associated with this exact subclonotype.
    #[prost(uint32, repeated, packed = "false", tag = "8")]
    pub umi_counts: ::prost::alloc::vec::Vec<u32>,
    /// Read counts of contigs associated with this exact subclonotype chain. The
    /// number of elements in this vector is equal to the number of barcodes
    /// associated with this exact subclonotype.
    #[prost(uint32, repeated, packed = "false", tag = "9")]
    pub read_counts: ::prost::alloc::vec::Vec<u32>,
    /// Names of contigs associated with this exact subclonotype chain. The number
    /// of elements in this vector is equal to the number of barcodes associated
    /// with this exact subclonotype. The contig name would be of the form
    /// `{barcode}_contig_{id}`.
    #[prost(string, repeated, tag = "10")]
    pub contig_ids: ::prost::alloc::vec::Vec<::prost::alloc::string::String>,
    /// Alignment of the `nt_sequence` to the nucleotide sequence of the clonotype
    /// consensus of this chain.
    /// TODO: Do we need amino acid alignment info?
    #[prost(message, required, tag = "11")]
    pub clonotype_consensus_aln: Alignment,
    /// Alignment of the `nt_sequence` to the nucleotide sequence of the
    /// concatenated donor reference of this chain (defined elsewhere in this
    /// file).
    /// TODO: Default donor reference to universal reference?
    #[prost(message, required, tag = "12")]
    pub donor_reference_aln: Alignment,
    /// Alignment of the `nt_sequence` to the nucleotide sequence of the
    /// concatenated universal reference of this chain (defined elsewhere in this
    /// file).
    #[prost(message, required, tag = "13")]
    pub universal_reference_aln: Alignment,
    /// Index of the start of the FWR1 sequence in the `nt_sequence`.
    #[prost(uint32, optional, tag = "14")]
    pub fwr1_start: ::core::option::Option<u32>,
    /// Index of the start of the CDR1 sequence in the `nt_sequence`.
    #[prost(uint32, optional, tag = "15")]
    pub cdr1_start: ::core::option::Option<u32>,
    /// Index of the start of the FWR2 sequence in the `nt_sequence`.
    #[prost(uint32, optional, tag = "16")]
    pub fwr2_start: ::core::option::Option<u32>,
    /// Index of the start of the CDR2 sequence in the `nt_sequence`.
    #[prost(uint32, optional, tag = "17")]
    pub cdr2_start: ::core::option::Option<u32>,
    /// Index of the start of the FWR3 sequence in the `nt_sequence`.
    #[prost(uint32, optional, tag = "18")]
    pub fwr3_start: ::core::option::Option<u32>,
    /// Index of the end of the FWR4 sequence in the `nt_sequence` (exclusive).
    #[prost(uint32, optional, tag = "19")]
    pub fwr4_end: ::core::option::Option<u32>,
    /// Nucleotide percent identity with the donor reference, outside junction region.
    #[prost(float, required, tag = "20")]
    pub dna_percent: f32,
    /// Amino acid percent identity with the donor reference, outside junction region.
    #[prost(float, required, tag = "21")]
    pub aa_percent: f32,
}
/// The chains in a clonotype are ordered an hence they have a unique index.
/// An exact subclonotype within a clonotype might not have all the chains that
/// are present in the clonotype. This structure stores the exact subclonotype
/// chain along with the `index` of the corresponding chain in the parent
/// clonotype.
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct ExactSubClonotypeChainInfo {
    /// The index of this chain in the parent clonotype
    #[prost(uint32, required, tag = "1")]
    pub index: u32,
    #[prost(message, required, tag = "2")]
    pub chain: ExactSubClonotypeChain,
}
/// Define an exact subclonotype.
///
/// All the barcodes within an exact subclonotype have the same number of
/// productive chains, the same sequence from the start of the V-REGION to the
/// end of the J-REGION as well as the same C-REGION annotation for each chain.
/// TODO: Maybe mutations outside V-J?
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct ExactSubClonotype {
    /// The chains in an exact subclonotype along with the index of the chain in
    /// the parent clonotype
    #[prost(message, repeated, tag = "1")]
    pub chains: ::prost::alloc::vec::Vec<ExactSubClonotypeChainInfo>,
    /// List of cell barcodes in this exact subclonotype. The number of elements in
    /// this list is equal to the number of elements in the `umi_counts` and
    /// `contig_ids` vector in an `ExactSubClonotypeChain`. The barcodes have the
    /// appropriate gem group as the suffix.
    #[prost(string, repeated, tag = "2")]
    pub cell_barcodes: ::prost::alloc::vec::Vec<::prost::alloc::string::String>,
    /// For a T cell, optionally annotate the exact subclonotype as one or both of
    /// iNKT/MAIT. This data structure stores the evidence we have about the
    /// annotation.  The evidence can be zero, and indeed this will be the case for
    /// all B cells.
    #[prost(message, required, tag = "3")]
    pub inkt_evidence: InvariantTCellAnnotation,
    #[prost(message, required, tag = "4")]
    pub mait_evidence: InvariantTCellAnnotation,
}
/// Define a clonotype chain
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct ClonotypeChain {
    /// The nucleotide sequence of this clonotype chain consensus
    /// What we actually compute here is not the consensus across the clonotype
    /// (whose biological meaning is questionable), but rather the sequence of the
    /// first exact subclonotype which has an entry for the given chain.  Over 99%
    /// of the time, this will be the first exact subclonotype.
    #[prost(bytes = "vec", required, tag = "1")]
    pub nt_sequence: ::prost::alloc::vec::Vec<u8>,
    /// Index of the 5' UTR region in the universal reference. The region in the
    /// universal reference is guaranteed to be `Region::U`
    #[prost(uint32, optional, tag = "2")]
    pub u_idx: ::core::option::Option<u32>,
    /// Index of the Variable region in the universal reference. The region in the
    /// universal reference is guaranteed to be `Region::V`
    #[prost(uint32, required, tag = "3")]
    pub v_idx: u32,
    /// Index of the Diversity region in the universal reference. The region in the
    /// universal reference is guaranteed to be `Region::D`. D-REGION is not
    /// present for light/alpha chains. Even for heavy/beta chains, this might be
    /// `None` if there is ambiguity in the annotation.
    #[prost(uint32, optional, tag = "4")]
    pub d_idx: ::core::option::Option<u32>,
    /// Index of the Joining region in the universal reference. The region in the
    /// universal reference is guaranteed to be `Region::J`
    #[prost(uint32, required, tag = "5")]
    pub j_idx: u32,
    /// Index of the Constant region in the universal reference. The region in the
    /// universal reference is guaranteed to be `Region::C`
    #[prost(uint32, optional, tag = "6")]
    pub c_idx: ::core::option::Option<u32>,
    /// Index of the Variable region in the donor reference. The region in the
    /// donor reference is guaranteed to be `Region::V` and the `universal_idx` in
    /// the donor reference item will be equal to the `v_idx`
    #[prost(uint32, optional, tag = "7")]
    pub donor_v_idx: ::core::option::Option<u32>,
    //// Index of the Joining region in the donor reference. The region in the
    //// donor reference is guaranteed to be `Region::J` and the `universal_idx` in
    //// the donor reference item will be equal to the `j_idx`
    #[prost(uint32, optional, tag = "8")]
    pub donor_j_idx: ::core::option::Option<u32>,
    /// Concatenated universal reference =
    ///     `nt_sequence` of universal_reference\[u_idx\] if u_idx is not None +
    ///     `nt_sequence` of universal_reference\[v_idx\] +
    ///     `nt_sequence` of universal_reference\[d_idx\] if d_idx is not None +
    ///     `nt_sequence` of universal_reference\[j_idx\] +
    ///     `nt_sequence` of universal_reference\[c_idx\] if c_idx is not None.
    #[prost(bytes = "vec", required, tag = "9")]
    pub universal_reference: ::prost::alloc::vec::Vec<u8>,
    /// Alignment of the `nt_sequence` to the nucleotide sequence of the
    /// concatenated universal reference of this chain.
    #[prost(message, required, tag = "10")]
    pub universal_reference_aln: Alignment,
    /// The concatenated donor reference is the same as the concatenated universal
    /// reference, however with substitutions:
    ///     `nt_sequence` of donor_reference\[donor_v_idx\] if donor_v_idx is not
    ///     None
    /// and
    ///     `nt_sequence` of donor_reference\[donor_j_idx\] if donor_j_idx is not
    ///     None.
    #[prost(bytes = "vec", required, tag = "11")]
    pub donor_reference: ::prost::alloc::vec::Vec<u8>,
    /// Alignment of the `nt_sequence` to the nucleotide sequence of the
    /// concatenated donor reference of this chain.
    #[prost(message, required, tag = "12")]
    pub donor_reference_aln: Alignment,
    /// Index of the start of the V-REGION in the `nt_sequence`.
    #[prost(uint32, required, tag = "13")]
    pub v_start: u32,
    /// Index of the stop of the aligned part of the V-REGION in the `nt-sequence`
    /// (exclusive).
    #[prost(uint32, required, tag = "14")]
    pub v_end: u32,
    /// Index of the stop of the aligned part of the V-REGION in universal
    /// reference V (exclusive).
    #[prost(uint32, required, tag = "15")]
    pub v_end_ref: u32,
    /// Index of the start of the aligned part of the J-REGION in the
    /// `nt_sequence`.
    #[prost(uint32, required, tag = "16")]
    pub j_start: u32,
    /// Index of the start of the aligned part of the J-REGION in universal
    /// reference J.
    #[prost(uint32, required, tag = "17")]
    pub j_start_ref: u32,
    /// Index of the end of the J-REGION in the `nt_sequence` (exclusive).
    #[prost(uint32, required, tag = "18")]
    pub j_end: u32,
    /// Index of the start of the CDR3 sequence in the `nt_sequence`. The start of
    /// the CDR3 amino acid in the `aa_sequence` is `(cdr3_start - v_start)/3`
    #[prost(uint32, required, tag = "19")]
    pub cdr3_start: u32,
    /// Index of the end of the CDR3 sequence in the `nt_sequence` (exclusive). The
    /// end of the CDR3 amino acid in the `aa_sequence` is `(cdr3_end - v_start)/3`
    #[prost(uint32, required, tag = "20")]
    pub cdr3_end: u32,
    /// IGH, IGK, IGL, TRA, TRB, TRD, TRG or conceivably something else.
    #[prost(string, required, tag = "21")]
    pub chain_type: ::prost::alloc::string::String,
    /// AA sequence of the clonotype chain consensus
    #[prost(bytes = "vec", required, tag = "22")]
    pub aa_sequence: ::prost::alloc::vec::Vec<u8>,
    /// AA sequence of the concatenated universal reference starting from the V
    /// regions
    #[prost(bytes = "vec", required, tag = "23")]
    pub aa_sequence_universal: ::prost::alloc::vec::Vec<u8>,
    /// AA sequence of the concatenated donor reference starting from the V regions
    #[prost(bytes = "vec", required, tag = "24")]
    pub aa_sequence_donor: ::prost::alloc::vec::Vec<u8>,
    /// Index of the start of the FWR1 sequence in the `nt_sequence`.
    #[prost(uint32, optional, tag = "25")]
    pub fwr1_start: ::core::option::Option<u32>,
    /// Index of the start of the CDR1 sequence in the `nt_sequence`.
    #[prost(uint32, optional, tag = "26")]
    pub cdr1_start: ::core::option::Option<u32>,
    /// Index of the start of the FWR2 sequence in the `nt_sequence`.
    #[prost(uint32, optional, tag = "27")]
    pub fwr2_start: ::core::option::Option<u32>,
    /// Index of the start of the CDR2 sequence in the `nt_sequence`.
    #[prost(uint32, optional, tag = "28")]
    pub cdr2_start: ::core::option::Option<u32>,
    /// Index of the start of the FWR3 sequence in the `nt_sequence`.
    #[prost(uint32, optional, tag = "29")]
    pub fwr3_start: ::core::option::Option<u32>,
    /// Index of the end of the FWR4 sequence in the `nt_sequence` (exclusive).
    #[prost(uint32, optional, tag = "30")]
    pub fwr4_end: ::core::option::Option<u32>,
}
/// Definition of a clonotype.
///
/// A clonotype is composed of a list of exact subclonotypes
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct Clonotype {
    /// The list of chains associated with this clonotype. The ordering of the
    /// chains is important as this order is preserved in the list of chains
    /// specified under each exact subclonotype. By convention heavy chain/beta
    /// chain comes ahead of light chain/alpha chain.
    /// TODO: What is the ordering when multiple chains of same kind are present?
    #[prost(message, repeated, tag = "1")]
    pub chains: ::prost::alloc::vec::Vec<ClonotypeChain>,
    /// The list of exact subclonotypes in this clonotype ordered by the number of
    /// cell barcodes in the exact subclonotype in descending order (TODO: Verify
    /// sort order). The number of chains listed under each exact subclonotype will
    /// be equal to the number of chains in this clonotype in the same order.
    /// However, some of the exact subclonotype chains could be `None`.
    #[prost(message, repeated, tag = "2")]
    pub exact_clonotypes: ::prost::alloc::vec::Vec<ExactSubClonotype>,
    /// The total number of cell barcodes associated with this clonotype. This can
    /// be inferred by summing up the number of barcodes within each exact
    /// subclonotype, but it is stored here for convenience.
    #[prost(uint32, required, tag = "3")]
    pub frequency: u32,
}
/// A single donor reference sequence and metadata packaged in a convenient
/// struct. In the current version of enclone, the donor reference is only
/// inferred for the V-REGIONs. But we could extend it to J-REGIONs in the
/// future.
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct UniversalReferenceItem {
    /// A unique identifier for this reference sequence that traces back to the
    /// reference fasta. Need not be ordered or continuous
    #[prost(uint32, required, tag = "1")]
    pub ref_idx: u32,
    /// The display name of this gene which will be shown in Loupe (optional allele
    /// information)
    #[prost(string, required, tag = "2")]
    pub display_name: ::prost::alloc::string::String,
    /// One of the U/V/D/J/C regions
    #[prost(enumeration = "Region", required, tag = "3")]
    pub region: i32,
    /// Nucleotide sequence associated with this reference item
    #[prost(bytes = "vec", required, tag = "4")]
    pub nt_sequence: ::prost::alloc::vec::Vec<u8>,
}
/// List of all universal reference sequences and metadata packaged in a
/// convenient struct. Currently just a list of items. Could add metadata like
/// species etc.
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct UniversalReference {
    /// UV(D)JC regions associated with a clonotype chain or an exact subclonotype
    /// chain are stored as indices into this vector
    #[prost(message, repeated, tag = "1")]
    pub items: ::prost::alloc::vec::Vec<UniversalReferenceItem>,
}
/// A single donor reference sequence and metadata packaged in a convenient
/// struct. In the current version of enclone, the donor reference is only
/// inferred for the V-REGIONs. But we could extend it to J-REGIONs in the
/// future.
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct DonorReferenceItem {
    /// Index of the parent sequence in the universal reference
    #[prost(uint32, required, tag = "1")]
    pub universal_idx: u32,
    /// Index of the donor associated with this reference. If there are no donors
    /// specified, this will be `0` by default.
    #[prost(uint32, required, tag = "2")]
    pub donor_idx: u32,
    /// The display name of this gene which will be shown in Loupe (optional allele
    /// information)
    /// TODO: Should this be modified to explicitly point out the donor? e.g TRAV-1
    /// [Donor 0]? for now, like this: "TRAV-1, donor 1, alt allele 1", etc.
    #[prost(string, required, tag = "3")]
    pub display_name: ::prost::alloc::string::String,
    /// Currently, the donor reference region will only be the V-REGION
    #[prost(enumeration = "Region", required, tag = "4")]
    pub region: i32,
    /// The nucleotide sequence associated with this reference item
    #[prost(bytes = "vec", required, tag = "5")]
    pub nt_sequence: ::prost::alloc::vec::Vec<u8>,
    /// Alignment of the `nt-sequence` with the nucleotide sequence of the
    /// corresponding universal reference item.
    #[prost(message, required, tag = "6")]
    pub universal_aln: Alignment,
}
/// List of all donor reference sequences and metadata packaged in a convenient
/// struct. Currently just a list of items. Could add metadata like all donor
/// names etc.
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct DonorReference {
    /// All the entries in this reference
    /// The donor V-REGION associated with a clonotype chain is stored as an index
    /// to this vector.
    #[prost(message, repeated, tag = "1")]
    pub items: ::prost::alloc::vec::Vec<DonorReferenceItem>,
}
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct GemWellInfo {
    #[prost(string, required, tag = "1")]
    pub donor: ::prost::alloc::string::String,
    #[prost(string, required, tag = "2")]
    pub origin: ::prost::alloc::string::String,
    #[prost(string, required, tag = "3")]
    pub library_id: ::prost::alloc::string::String,
    /// Data for all the additional columns. For convenience, we are storing it as
    /// a map, rather than an array with the order defined by `additional_columns`
    /// in Metadata
    #[prost(map = "string, string", tag = "8")]
    pub additional_data:
        ::std::collections::HashMap<::prost::alloc::string::String, ::prost::alloc::string::String>,
}
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct Metadata {
    /// The additional metadata columns specified by the user
    #[prost(string, repeated, tag = "1")]
    pub additional_columns: ::prost::alloc::vec::Vec<::prost::alloc::string::String>,
    /// List of all the donors. The `donor_idx` in the `DonorReference` is an index
    /// into this array. This will be empty for a single sample case
    #[prost(string, repeated, tag = "2")]
    pub donors: ::prost::alloc::vec::Vec<::prost::alloc::string::String>,
    /// Key: Gem well
    /// Value: Metadata for each gem well
    /// This will be empty for a single sample case
    #[prost(map = "uint32, message", tag = "3")]
    pub per_gem_well_info: ::std::collections::HashMap<u32, GemWellInfo>,
}
/// Outputs from a single enclone run.
///
/// This message itself is not written in the proto file, but the order of
/// messages follow the order of fields in this message
#[derive(::serde::Serialize, ::serde::Deserialize, Clone, PartialEq, ::prost::Message)]
pub struct EncloneOutputs {
    #[prost(string, required, tag = "1")]
    pub version: ::prost::alloc::string::String,
    /// Metadata associated with this run. For single sample runs, there is no
    /// useful metadata. In the context of aggr, this will be populated from the
    /// user input
    #[prost(message, required, tag = "10")]
    pub metadata: Metadata,
    /// The universal reference used.
    #[prost(message, required, tag = "20")]
    pub universal_reference: UniversalReference,
    /// The inferred donor reference.
    #[prost(message, required, tag = "30")]
    pub donor_reference: DonorReference,
    /// The total number of clonotypes
    #[prost(uint32, required, tag = "100")]
    pub num_clonotypes: u32,
    /// List of all clonotypes computed in this enclone run. Each clonotype is
    /// stored as an individual message in order to enable streaming.
    #[prost(message, repeated, tag = "110")]
    pub clonotypes: ::prost::alloc::vec::Vec<Clonotype>,
}
/// Various regions within a VDJ transcript
#[derive(
    ::serde::Serialize,
    ::serde::Deserialize,
    Clone,
    Copy,
    Debug,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
    ::prost::Enumeration,
)]
#[repr(i32)]
pub enum Region {
    /// 5' untranslated region
    U = 0,
    /// Variable region
    V = 1,
    /// Diversity region
    D = 2,
    /// Joining region
    J = 3,
    /// Constant region
    C = 4,
}
