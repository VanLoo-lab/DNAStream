def test_add_snv_trees(dnastream_obj, tree_file):
    """Test adding SNV trees and safe mode behavior."""

    # Add trees for the first time
    dnastream_obj.add_trees_from_file(tree_file, tree_type="SNV", method="conipher")
    tree_data = dnastream_obj._get_data("tree/SNV_trees/data")

    assert tree_data.shape[0] > 0, "No trees were added!"
    initial_tree_count = tree_data.shape[0]

    # Ensure safe mode prevents duplicate additions
    dnastream_obj.add_trees_from_file(tree_file, tree_type="SNV", method="conipher")
    tree_data = dnastream_obj._get_data("tree/SNV_trees/data")

    assert (
        tree_data.shape[0] == initial_tree_count
    ), "Safe mode failed to prevent duplicate trees!"

    # Disable safe mode and append trees
    dnastream_obj.safe_mode_disable()
    dnastream_obj.add_trees_from_file(tree_file, tree_type="SNV", method="conipher")
    tree_data = dnastream_obj._get_data("tree/SNV_trees/data")

    assert (
        tree_data.shape[0] == initial_tree_count * 2
    ), "Disabling safe mode did not allow appending new trees!"


def test_add_snv_trees_from_edge_list(dnastream_obj):
    """Test adding SNV trees from an edge list."""
    edge_list = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]
    tree_list = [edge_list]
    dnastream_obj.add_trees_from_edge_lists(tree_list, method="conipher")
    tree_data = dnastream_obj._get_data("tree/SNV_trees/data")

    assert (
        tree_data.shape[0] == len(tree_list)
    ), "Failed to add trees from edge lists!"
  

    # # test sapling
    # dna_stream_obj.safe_mode_enable()
    # dna_stream_obj.add_trees_from_file(SAPLING, tree_type="SNV", method="sapling")
    # tree_dat = dna_stream_obj._get_data("tree/SNV_trees/data")
    # assert tree_dat.shape[0] > numtrees * 2
