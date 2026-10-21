pub fn runs(labels: &[u32]) -> Vec<(u32, u32)> {
    let mut sorted = labels.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    let mut runs: Vec<(u32, u32)> = Vec::new();
    for label in sorted {
        match runs.last_mut() {
            Some((_, end)) if end.checked_add(1) == Some(label) => *end = label,
            _ => runs.push((label, label)),
        }
    }
    runs
}

pub fn write(labels: &[u32]) -> String {
    runs(labels)
        .into_iter()
        .map(|(start, end)| match start == end {
            true => start.to_string(),
            false => format!("{start}..{end}"),
        })
        .collect::<Vec<String>>()
        .join(",")
}

pub fn span(chain: &str, start: u32, end: u32) -> String {
    match start == end {
        true => format!("{chain}{start}"),
        false => format!("{chain}{start}-{end}"),
    }
}

pub fn contig(chain: &str, labels: &[u32]) -> String {
    runs(labels)
        .into_iter()
        .map(|(start, end)| span(chain, start, end))
        .collect::<Vec<String>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn consecutive_labels_collapse_into_one_run() {
        assert_eq!(runs(&[1, 2, 3, 4]), [(1, 4)]);
        assert_eq!(write(&[1, 2, 3, 4]), "1..4");
        assert_eq!(contig("A", &[1, 2, 3, 4]), "A1-4");
    }

    #[test]
    fn gaps_start_a_new_run() {
        assert_eq!(runs(&[1, 2, 5, 6, 9]), [(1, 2), (5, 6), (9, 9)]);
        assert_eq!(write(&[1, 2, 5, 6, 9]), "1..2,5..6,9");
        assert_eq!(contig("B", &[1, 2, 5, 6, 9]), "B1-2,B5-6,B9");
    }

    #[test]
    fn a_span_is_written_from_its_ends() {
        assert_eq!(span("A", 1, 3), "A1-3");
        assert_eq!(span("A", 4, 4), "A4");
    }

    #[test]
    fn a_single_label_needs_no_range() {
        assert_eq!(write(&[7]), "7");
        assert_eq!(contig("A", &[7]), "A7");
    }

    #[test]
    fn repeated_labels_do_not_extend_a_run() {
        assert_eq!(runs(&[3, 3, 4]), [(3, 4)]);
        assert_eq!(write(&[]), "");
    }

    #[test]
    fn unsorted_labels_are_ordered_before_they_are_written() {
        assert_eq!(write(&[5, 6, 1, 2]), "1..2,5..6");
        assert_eq!(write(&[2, 1, 1]), "1..2");
        assert_eq!(contig("A", &[9, 3, 4]), "A3-4,A9");
    }

    #[test]
    fn the_highest_label_does_not_overflow_a_run() {
        assert_eq!(runs(&[u32::MAX, 4]), [(4, 4), (u32::MAX, u32::MAX)]);
        assert_eq!(runs(&[u32::MAX - 1, u32::MAX]), [(u32::MAX - 1, u32::MAX)]);
        assert_eq!(write(&[u32::MAX]), "4294967295");
    }
}
