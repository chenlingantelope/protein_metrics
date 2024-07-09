ESM_IF = True  # @param {type:"boolean"}
ProteinMPNN = True  # @param {type:"boolean"}
MIF_ST = True  # @param {type:"boolean"}
AlphaFold2_pLDDT = True  # @param {type:"boolean"}

# ESM-IF
import esm
import esm.inverse_folding
from utils import map_coords, aligned_to_raw_seq
esm_if_model, esm_if_alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
esm_if_model.eval()

def esm_if_score(pdb_seq, design_seq, coords):
    raw_pdb_seq = aligned_to_raw_seq(str(pdb_seq.seq))
    raw_design_seq = aligned_to_raw_seq(str(design_seq.seq))
    mapped_coords = map_coords(coords, design_seq.seq, pdb_seq.seq)
    _, ll_pdb = esm.inverse_folding.util.score_sequence(
        esm_if_model, esm_if_alphabet, coords, raw_pdb_seq)
    _, ll_design = esm.inverse_folding.util.score_sequence(
        esm_if_model, esm_if_alphabet, mapped_coords, raw_design_seq)
    return ll_pdb, ll_design


# the following can be run in main.py to visualized the loss across the sequence
# loss_pdb = esm.inverse_folding.util.get_sequence_loss(
#     esm_if_model, esm_if_alphabet, coords, raw_pdb_seq)
#
# loss_design = esm.inverse_folding.util.get_sequence_loss(
#     esm_if_model, esm_if_alphabet, mapped_coords, raw_design_seq)
#
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots(figsize=(20, 8))
# ax.plot(loss_pdb[0], label="Loss")
# ax.plot(loss_design[0], label="Design Loss")
# ax.legend()
# plt.show()


# # ProteinMPNN
# if ProteinMPNN:
#     with tempfile.TemporaryDirectory() as output_dir:
#         for i, pdb_file in enumerate(glob("/content/pdbs/*.pdb")):
#             command_line_arguments = [
#                 "python",
#                 "ProteinMPNN/vanilla_proteinmpnn/protein_mpnn_run.py",
#                 "--pdb_path", pdb_file,
#                 "--pdb_path_chains", "A",
#                 "--score_only", "1",
#                 "--save_score", "1",
#                 "--out_folder", output_dir,
#                 "--batch_size", "1"
#             ]
#             fstem = Path(pdb_file).stem
#             name = fstem
#             outfile = output_dir + f"outfile_{i}.txt"
#             with open(outfile, "w") as fh:
#                 proc = subprocess.run(command_line_arguments, stdout=subprocess.PIPE, check=True)
#                 print(proc.stdout.decode('utf-8'), file=fh)
#             with open(outfile, "r") as score_file_h:
#                 score_file_lines = score_file_h.readlines()
#             score_line = score_file_lines[-2].split(",")
#             score_parts = score_line[1].strip().split(": ")
#             assert score_parts[0] == "mean"
#             score = -1 * float(score_parts[1])
#             add_metric(results, name, "ProteinMPNN", score)
#
# # MIF-ST
#
# if MIF_ST:
#     with tempfile.TemporaryDirectory() as output_dir:
#         spec_file_path = output_dir + "/spec_file.tsv"
#         with open(spec_file_path, 'w') as f:
#             f.write('name\tsequence\tpdb\n')
#             for pdb_file in glob("/content/pdbs/*.pdb"):
#                 seq = get_pdb_sequence(pdb_file)
#                 name = Path(pdb_file).stem
#                 f.write(name + '\t' + seq + '\t' + pdb_file + '\n')
#         # print(spec_file_path)
#         proc = subprocess.run(
#             ['python', "/tmp/extract_mif.py", "mifst", spec_file_path, output_dir + "/", "logits", "--include", "logp",
#              "--device", device], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
#         # print(proc.stderr.decode("utf-8"))
#         # print(proc.stdout.decode("utf-8"))
#         df = pd.read_table(output_dir + '/mifst_logp.tsv')
#         df = df.rename(columns={'name': 'id', 'logp': 'mifst_logp'}, )
#         for _, row in df.iterrows():
#             add_metric(results, row["id"], "MIF-ST", row["mifst_logp"])
#
# # pLDDT
# if AlphaFold2_pLDDT:
#     for pdb_file in glob("/content/pdbs/*.pdb"):
#         fstem = Path(pdb_file).stem
#         name = fstem
#         pdb_file = pdb.PDBFile.read(pdb_file)
#         atoms = pdb_file.get_structure(extra_fields=['b_factor'])
#         prev_residue = -1
#         plddt_sum = 0
#         residue_count = 0
#         for a in atoms[0]:
#             if a.res_id != prev_residue:
#                 prev_residue = a.res_id
#                 residue_count += 1
#                 plddt_sum += a.b_factor
#         add_metric(results, name, "AlphaFold2 pLDDT", plddt_sum / residue_count)
