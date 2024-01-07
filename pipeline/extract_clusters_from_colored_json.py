import json
import os
from tqdm import tqdm

if __name__ == "__main__":
    file_path = r"/Users/luca/dev/GPCRs for HMMBuild/output/panphylum_correct18052022/Correct/rooted_with_iTOL/cluster_TF315321_colored.json"
    target_folder = file_path.replace(".json", "") + "_clusters"
    os.makedirs(target_folder, exist_ok=True)
    with open(file_path, "r") as fp:
        out = json.load(fp)
        node_names = out.keys()
        for n in tqdm(node_names, desc="Dumping entries per colored cluster"):
            base_name = os.path.join(target_folder, 
                                     os.path.basename(file_path).replace("colored", n))
            ext = out[n]["cluster"]
            if ext:
                with open(base_name, "w") as fw:
                    json.dump(ext, fw)
            else:
                print(f"No cluster found for {n}.")
