def upgma(distances):
    clusters = {name: [name] for name in distances.keys()}
    heights = {name: 0 for name in distances.keys()}

    step = 1

    while len(clusters) > 1:
        # Find closest pair of clusters
        min_dist = float("inf")
        pair = None

        names = list(clusters.keys())
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = names[i], names[j]
                d = distances[a][b]

                if d < min_dist:
                    min_dist = d
                    pair = (a, b)

        a, b = pair
        new_cluster = f"({a},{b})"

        print(f"\nStep {step}: merge {a} and {b}")
        print(f"Distance = {min_dist}")
        print(f"Branch height = {min_dist / 2}")

        # Create new cluster
        clusters[new_cluster] = clusters[a] + clusters[b]
        heights[new_cluster] = min_dist / 2

        # Add distances for new cluster
        distances[new_cluster] = {}

        for other in list(clusters.keys()):
            if other not in [a, b, new_cluster]:
                size_a = len(clusters[a])
                size_b = len(clusters[b])

                new_dist = (
                    distances[a][other] * size_a +
                    distances[b][other] * size_b
                ) / (size_a + size_b)

                distances[new_cluster][other] = new_dist
                distances[other][new_cluster] = new_dist

        distances[new_cluster][new_cluster] = 0

        # Remove old clusters
        del clusters[a]
        del clusters[b]

        for key in list(distances.keys()):
            distances[key].pop(a, None)
            distances[key].pop(b, None)

        distances.pop(a, None)
        distances.pop(b, None)

        print("Current clusters:", list(clusters.keys()))
        step += 1

    return list(clusters.keys())[0]


distances = {
    "S1": {"S1": 0, "S2": 3, "S3": 4, "S4": 5},
    "S2": {"S1": 3, "S2": 0, "S3": 6, "S4": 8},
    "S3": {"S1": 4, "S2": 6, "S3": 0, "S4": 9},
    "S4": {"S1": 5, "S2": 8, "S3": 9, "S4": 0},
}

tree = upgma(distances)

print("\nFinal UPGMA tree:")
print(tree)

