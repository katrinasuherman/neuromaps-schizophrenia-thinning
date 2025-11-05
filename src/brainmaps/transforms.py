from neuromaps.transforms import fsaverage_to_fslr, fslr_to_fslr, civet_to_fslr

def to_fsLR32k_from_fsaverage10k(annotation):
    return fsaverage_to_fslr(annotation, target_density="32k")

def to_fsLR32k_from_fsLR164k(annotation, hemi=None):
    return fslr_to_fslr(annotation, target_density="32k", hemi=hemi)

def to_fsLR32k_from_civet41k(annotation):
    return civet_to_fslr(annotation, target_density="32k")
