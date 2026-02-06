#!/bin/python3
import json
import os
import sys

if len(sys.argv) < 2:
    print("Usage: python script.py <sample_path>")
    sys.exit(1)

sample = sys.argv[1]
json_file = f"{sample}"

if os.path.exists(json_file):
    with open(json_file) as f:
        data = json.load(f)
    
    records = data.get('records', [])
    total_clusters = sum(len(r.get('areas', [])) for r in records)
    
    print(f"Total BGCs detected: {total_clusters}\n")
    
    # Count by type
    cluster_types = {}
    for record in records:
        for area in record.get('areas', []):
            for product in area.get('products', []):
                cluster_types[product] = cluster_types.get(product, 0) + 1
    
    if cluster_types:
        print("BGC types:")
        for ctype, count in sorted(cluster_types.items()):
            print(f"  {ctype}: {count}")
else:
    print(f"antiSMASH JSON not found: {json_file}")