(defun emacs-format-function ()
  "Format the whole buffer."
  (if (< 1 (count-windows))
      (delete-other-windows (selected-window)))
  (catch 'tag
    (while t
      (c-mode)
      (setq indent-tabs-mode 'nil)
      ;;(setq c-basic-offset 4)
      ;; can add a save hook to delete trailing whitespaces
      (delete-trailing-whitespace)
      (untabify (point-min) (point-max))
      (indent-region (point-min) (point-max) nil)
      (if buffer-file-name  ; nil for *scratch* buffer
          (progn
            (write-file buffer-file-name)
            (kill-buffer (current-buffer)))
        (throw 'tag t)))))
