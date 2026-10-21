const INDICATORS: [char; 20] = [
    '-', '?', ':', ',', '[', ']', '{', '}', '#', '&', '*', '!', '|', '>', '\'', '"', '%', '@', '`',
    '~',
];

const FLOW: [char; 5] = [',', '[', ']', '{', '}'];

pub fn yaml(text: &str) -> String {
    match plain(text) {
        true => text.to_string(),
        false => escaped(text, false),
    }
}

pub fn yaml_flow(text: &str) -> String {
    match plain(text) && !text.contains(FLOW) {
        true => text.to_string(),
        false => escaped(text, false),
    }
}

pub fn json(text: &str) -> String {
    escaped(text, true)
}

fn plain(text: &str) -> bool {
    !text.is_empty()
        && text.trim() == text
        && !text.starts_with(INDICATORS)
        && !text.contains(": ")
        && !text.contains(" #")
        && !text.ends_with(':')
        && !text.chars().any(char::is_control)
        && !retyped(text)
}

fn retyped(text: &str) -> bool {
    matches!(
        text.to_ascii_lowercase().as_str(),
        "null" | "true" | "false" | "yes" | "no" | "on" | "off"
    ) || text.parse::<f64>().is_ok()
        || text.parse::<i64>().is_ok()
}

fn escaped(text: &str, json: bool) -> String {
    let mut out = String::with_capacity(text.len() + 2);
    out.push('"');
    for letter in text.chars() {
        match letter {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            letter if letter.is_control() => {
                let code = letter as u32;
                match (json, code < 0x100) {
                    (true, _) => out.push_str(&format!("\\u{code:04x}")),
                    (false, true) => out.push_str(&format!("\\x{code:02x}")),
                    (false, false) => out.push_str(&format!("\\u{code:04x}")),
                }
            }
            letter => out.push(letter),
        }
    }
    out.push('"');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ordinary_names_stay_plain() {
        for text in [
            "A",
            "S",
            "Q",
            "WHL",
            "SG",
            "target.cif",
            "10..14C6C3",
            "1..6",
            "4..5",
        ] {
            assert_eq!(yaml(text), text, "{text} should stay plain");
            assert_eq!(yaml_flow(text), text, "{text} should stay plain in flow");
        }
    }

    #[test]
    fn indicator_characters_are_quoted() {
        assert_eq!(yaml("*"), "\"*\"");
        assert_eq!(yaml("[odd].cif"), "\"[odd].cif\"");
        assert_eq!(yaml("{a}.cif"), "\"{a}.cif\"");
        assert_eq!(yaml(",leading.cif"), "\",leading.cif\"");
        assert_eq!(yaml("-dash"), "\"-dash\"");
        assert_eq!(yaml("@at"), "\"@at\"");
    }

    #[test]
    fn values_that_would_change_type_are_quoted() {
        for text in [
            "no", "No", "NO", "true", "off", "null", "~", "1", "1.5", "1e5", "0",
        ] {
            let written = yaml(text);
            assert!(
                written.starts_with('"'),
                "{text} should be quoted, got {written}"
            );
        }
        assert_eq!(yaml("n"), "n");
        assert_eq!(yaml("y"), "y");
    }

    #[test]
    fn whitespace_and_comments_are_quoted() {
        assert_eq!(yaml(" leading"), "\" leading\"");
        assert_eq!(yaml("trailing "), "\"trailing \"");
        assert_eq!(yaml("a: b"), "\"a: b\"");
        assert_eq!(yaml("a #b"), "\"a #b\"");
        assert_eq!(yaml("ends:"), "\"ends:\"");
        assert_eq!(yaml(""), "\"\"");
    }

    #[test]
    fn control_characters_are_escaped_for_each_format() {
        assert_eq!(yaml("a\tb"), "\"a\\tb\"");
        assert_eq!(yaml("a\nb"), "\"a\\nb\"");
        assert_eq!(yaml("a\u{1}b"), "\"a\\x01b\"");
        assert_eq!(json("a\tb"), "\"a\\tb\"");
        assert_eq!(json("a\u{1}b"), "\"a\\u0001b\"");
        assert_eq!(json("a\"b\\c"), "\"a\\\"b\\\\c\"");
    }

    #[test]
    fn flow_values_quote_their_separators() {
        assert_eq!(yaml_flow("a,b"), "\"a,b\"");
        assert_eq!(yaml_flow("a]b"), "\"a]b\"");
        assert_eq!(yaml("a,b"), "a,b");
    }
}
