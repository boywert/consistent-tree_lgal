#ifndef FIND_PARENTS_AND_CLEANUP_H
#define FIND_PARENTS_AND_CLEANUP_H

void cleanup_print_halos(int64_t output_num, int64_t stage);
void cleanup_phantoms(float a1, float a2, int stage);
inline int64_t id_to_index(struct halo_stash h, int64_t id);
void find_parents(void);
int sort_halo_order(const void *a, const void *b);
void clear_halo_stash(struct halo_stash *h);
void fix_halos_without_descendants(float a1, float a2);
void zero_halo_stash(struct halo_stash *h);
void cleanup_find_new_descendants(float a1, float a2);
void cleanup_short_tracks(float a1, float a2, int stage, int64_t output_index, int64_t num_outputs);
void translate_ids(int64_t stage);
void detect_major_mergers(float a);
void tag_children_of_centrals(float a);


#endif /* FIND_PARENTS_AND_CLEANUP_H */
